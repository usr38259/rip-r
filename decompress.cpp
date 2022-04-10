#include "decompress.h"
#include <stack>

using namespace std;

vector<byte> decompress(const vector<byte>& data)
{
	auto inptr = data.cbegin();
	int bitbuf = 1;
	auto readbit = [&]() {
		if (bitbuf == 1) {
			if (inptr == data.end()) {
				throw exception("Unexpected end of data");
			}
			bitbuf = (*inptr++) + 0x100;
		}
		int cf = bitbuf & 1;
		bitbuf >>= 1;
		return cf;
	};

	auto readFromTree = [&](vector<int>& tree) {
		int i = 0;
		while (true) {
			if (readbit() == 1) {
				i++;
			}
			if (i >= (int)tree.size()) {
				throw exception("Invalid compressed data");
			}
			i = tree[i];
			if (i < 0) {
				return ~i;
			}
		}
	};

	auto buildTree = [&](vector<int>& depths, int offset, int count) {
		int q = 0;
		for(int i=0; i < count; i++) {
			if (depths[offset+i] > 0) q++;
		}
		// Building tree of q leafs.
		// Root node is omitted, so the total number of nodes must be q*2-2.
		if (q < 2) {
			throw exception("Invalid compressed data");
		}
		vector<int> tree;
		int childIdx = 0;
		int nodeIdx = 0, nodeDepth = 1;
		stack<int> st;
		st.push(0); // root node depth
		st.push(nodeIdx);
		st.push(nodeDepth);
		for (int curDepth = 1; ; curDepth++)
		{
			if (curDepth > 15 || (int)tree.size() > q*2-2) {
				throw exception("Invalid compressed data");
			}
			for (int code = count; code > 0; ) {
				do {
					code--;
				} while (code >= 0 && depths[offset + code] != curDepth);
				if (code < 0) {
					break;
				}
				while (nodeDepth < curDepth) {
					childIdx += 2;
					while((int)tree.size() < nodeIdx+1) tree.push_back(0);
					tree[nodeIdx] = childIdx; // ptr to children
					nodeIdx = childIdx;
					nodeDepth++;
					st.push(nodeIdx);
					st.push(nodeDepth);
				}
				while((int)tree.size() < nodeIdx+1) tree.push_back(0);
				tree[nodeIdx] = ~code; // leaf
				nodeDepth = st.top(); st.pop();
				if (nodeDepth == 0) {
					break;
				}
				nodeIdx = st.top(); st.pop();
				nodeIdx++;
			}
			if (nodeDepth == 0) {
				break;
			}
		}

		if (tree.size() != q*2-2) {
			throw exception("Invalid compressed data");
		}
		return tree;
	};

	vector<int> depths1(18);
	for (int i = 0; i < 18; i++) {
		depths1[i] = readbit()*8 + readbit()*4 + readbit()*2 + readbit();
	}
	auto tree1 = buildTree(depths1, 0, 18);

	vector<int> depths2; // may contain up to 0x142 items, but only 0x140 used
	while (true) {
		if (depths2.size() > 0x120 + 32) {
			throw exception("Invalid compressed data");
		}
		int code = readFromTree(tree1);
		if (code < 16) {
			depths2.push_back(code);
		} else { // 16, 17
			if (depths2.empty()) {
				throw exception("Invalid compressed data");
			}
			int t = depths2.back();
			depths2.push_back(t);
			depths2.push_back(t);
			if (code == 16) {
				break;
			}
		}
	}

	if (depths2.size() < 0x120 + 32) {
		throw exception("Invalid compressed data");
	}

	vector<int> codesAndCountTree = buildTree(depths2, 0, 0x120);
	vector<int> offsetTree = buildTree(depths2, 0x120, 32);

	vector<byte> result;

	int offset = 0;
	while (true)
	{
		int code = readFromTree(codesAndCountTree);
		if (code < 0x100) {
			result.push_back((byte)code);
		}
		else if (code == 0x100) {
			break;
		}
		else {
			code &= 0xFF;
			int count = code; // 01..1f
			if (code >= 5) {
				code = (code + 0xFB + 2) & 0xFF; // 02..1c
				count = 2 + (code & 1); // 2..3
				int extrabits = code / 2; // 01..0e
				for (int i = 0; i < extrabits; i++) {
					count = count * 2 + readbit();
				}
				count++; // 5..49152
			}

			code = readFromTree(offsetTree); // 00..1f
			if (code != 0) {
				offset = code;  // 01..1f
				if (code >= 5) {
					code = (code + 0xFB + 2) & 0xFF; // 02..1c
					offset = 2 + (code & 1); // 2..3
					int extrabits = code / 2; // 01..0e
					for (int i = 0; i < extrabits; i++) {
						offset = offset * 2 + readbit();
					}
					offset++; // 5..49152
				}

				if (offset >= 0x100) {
					count++;
				}
			}

			if (count < 1 || offset < 1 || (size_t)offset > result.size()) {
				throw exception("Invalid compressed data");
			}
			for (int i = 0; i < count; i++) {
				byte t = result[result.size() - offset];
				result.push_back(t);
			}
		}
	}

	if (inptr != data.end()) {
		throw exception("Extra data in the end");
	}

	return result;
};

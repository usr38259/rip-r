/*
 * Copyright © 2022 Eugene Larchenko <el6345@gmail.com>. All rights reserved.
 * See the attached LICENSE file.
 */

#include "compress.h"
#include "types.h"
#include <stdio.h> // printf
#include <vector>
#include <string>
#include <stack>
#include <functional>
#include <time.h> // clock
//#include <algorithm> // min,max
#include "longest_match.h"

using namespace std;

const int MAX_CODE_LEN = 15;
const int EOS_CODE = 256;
const int MAX_OFFSET = 49152;
const int MAX_LEN = 49152; // 49153 for offset>=256
const int MAX_LEN_CORRECTION = 1;

// ranges for codes 0..31
const ushort range_floor[32] = { 0, 1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 0x101, 0x181, 0x201, 0x301, 0x401, 0x601, 0x801, 0xc01, 0x1001, 0x1801, 0x2001, 0x3001, 0x4001, 0x6001, 0x8001 };
const ushort range_ceil[32] = { 0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 0x180, 0x200, 0x300, 0x400, 0x600, 0x800, 0xc00, 0x1000, 0x1800, 0x2000, 0x3000, 0x4000, 0x6000, 0x8000, 0xc000 };

double timeSpent;

void ZFunction(const byte* zdata, int startpos, int len, int* matchLen)
{
	if (startpos >= len) {
		throw exception();
	}

	int* z = matchLen;

	const byte* s = zdata;
	const int sstart = startpos;
	const int n = (len - startpos) + len;

	int l = 0, r = 0;
	for (int i = 1; i < len; i++)
	{
		int zi;
		if (i > r)
			zi = 0;
		else
			zi = min(z[i - l], r + 1 - i);

		while (i + zi < n && s[sstart + i + zi] == s[sstart + zi]) {
			zi++;
		}

		if (i + zi - 1 > r)	{
			r = i + zi - 1;
			l = i;
		}

		z[i] = zi;
	}

	// z[0] is undefined
}

struct BitStream
{
	size_t Count;
	vector<byte> Data;

	BitStream()	{
		Count = 0;
		Data.clear();
	}

	void Write(int count, int bits) {
		if ((bits >> count) != 0) {
			throw exception();
		}
		for (int i = count-1; i >= 0; i--) {
			if (this->Count % 8 == 0) {
				Data.push_back(0);
			}
			int v = bits >> i & 1;
			Data[this->Count / 8] |= (v << this->Count % 8);
			this->Count++;
		}
	}

	void Write(string bits) {
		for (char v : bits) {
			if (this->Count % 8 == 0) {
				Data.push_back(0);
			}
			Data[this->Count / 8] |= ((v & 1) << this->Count % 8);
			this->Count++;
		}
	}

	int operator[](size_t i) {
		return Data[i / 8] >> (i % 8) & 1;
	}
};

struct rip_tree_builder
{
	static vector<int> BuildTable(const vector<int>& depths, int offset, int count)
	{
		int q = 0;
		for(int i=0; i < count; i++) {
			if (depths[offset+i] > 0) q++;
		}
		// Building tree of q leafs.
		// Root node is omitted, so the total number of nodes must be q*2-2
		if (q < 2) {
			throw exception("too few codes");
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
			if (curDepth > MAX_CODE_LEN || (int)tree.size() > q*2-2) {
				throw exception("Invalid depths table");
			}
			for (int code = count; code > 0; )
			{
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
			throw exception("Invalid depths table");
		}
		return tree;
	}

	static vector<string> MakeCodes(const vector<int>& depths)
	{
		vector<int> tree = BuildTable(depths, 0, (int)depths.size());
		vector<string> res(depths.size());
		std::function<void(int,string)> traverse = [&](int n, string& code) {
			if (tree[n] < 0) {
				res[~tree[n]] = code;
			}
			else {
				traverse(tree[n] + 0, code + "0");
				traverse(tree[n] + 1, code + "1");
			}
		};
		// no root node, so traverse both its children
		traverse(0, "0");
		traverse(1, "1");
		return res;
	}
};


struct Node;
struct Node
{
	int Chr, Count;
	Node *L, *R;
	bool IsLeaf;

	Node() {
		memset(this, 0, sizeof(*this));
	}

	Node(int chr, int count) {
		Chr = chr; Count = count;
		IsLeaf = true;
	}

	Node(Node* l, Node* r, int count) {
		Chr = -1;
		Count = count;
		L = l; R = r;
		IsLeaf = false;
	}

	static int64 Compare_byCountAsc_thenByCodeDesc(Node& n1, Node& n2) {
		int64 t = (int64)n1.Count - n2.Count;
		if (t == 0) {
			if (n1.IsLeaf && n2.IsLeaf)
				t = -(n1.Chr - n2.Chr); // smaller codepoints tend to get shorter codes
		}
		return t;
	}

	static int64 Compare_byCountDesc_thenByCodeAsc(Node& n1, Node& n2) {
		int64 t = -((int64)n1.Count - n2.Count);
		if (t == 0) {
			if (n1.IsLeaf && n2.IsLeaf)
				t = (n1.Chr - n2.Chr); // smaller codepoints tend to get shorter codes
		}
		return t;
	}

};

vector<int> MakeCodeTable_LimitedLen(const vector<int>& stats, int codeLenLimit)
{
	// сортировка по count обязательна, без неё результат существенно хуже

	vector<Node> a;
	for(int c=0; c < (int)stats.size(); c++) {
		if (stats[c] > 0) {
			Node n(c, stats[c]);
			size_t i = a.size(); // insertion sort
			while(i>0 && Node::Compare_byCountDesc_thenByCodeAsc(a[i-1], n) > 0) {
				i--;
			}
			a.insert(a.begin()+i, n);
		}
	}

	if (a.size() == 0)
		throw exception("a.Length==0");

	const float infinity = 1e38f;
	const int max_a_size = 0x120;
	int a_size = (int)a.size();
	if (a_size > max_a_size || codeLenLimit > MAX_CODE_LEN) {
		throw exception();
	}

	auto dp = new float[max_a_size+1][max_a_size][MAX_CODE_LEN+2];
	auto dp2 = new short[max_a_size+1][max_a_size][MAX_CODE_LEN+2];
	if (!dp || !dp2) {
		throw bad_alloc();
	}
	for (int range = 1; range <= a_size; range++) {
		for (int start = 0; start <= a_size - range; start++)	{
			for (int depth = 0; depth <= codeLenLimit+1; depth++) {
				int bestm = -1;
				float res;
				if (depth > codeLenLimit) {
					res = infinity;
				}
				else if (range == 1) {
					res = (float)depth * a[start].Count;
				}
				else {
					res = infinity;
					for (int m = 1; m < range; m++) {
						// try split at position m
						auto t = dp[m][start][depth + 1] + dp[range - m][start + m][depth + 1];
						if (t < res) {
							res = t; bestm = m;
						}
					}
				}
				dp2[range][start][depth] = (short)bestm;
				dp[range][start][depth] = res;
			}
		}
	}
	auto reslen = dp[a_size][0][0]; // encoded size in bits
	if (reslen >= infinity)
		throw exception(); // no solution; can't happen

	vector<int> code_lengths(stats.size());
	std::function<void(int,int,int)> traverse = [&](int left, int range, int depth) {
		if (depth > codeLenLimit) {
			throw exception();
		}
		if (range == 1) {
			code_lengths[a[left].Chr] = depth;
		}
		else {
			int m = dp2[range][left][depth];
			traverse(left, m, depth + 1);
			traverse(left + m, range - m, depth + 1);
		}
	};
	traverse(0, a_size, 0);

	delete[] dp;
	delete[] dp2;

	return code_lengths;
}

vector<string> MakeCodeTable(vector<int> stats, int codeLenLimit)
{
	vector<Node> nodes(stats.size() * 2); int nq;
	vector<Node*> a;

	reconstruct_a:
	a.clear();
	nq = 0;
	for(int c=0; c < (int)stats.size(); c++) {
		if (stats[c] > 0) {
			nodes[nq++] = Node(c, stats[c]);
			Node* n = &nodes[nq-1];
			size_t i = a.size(); // insertion sort
			while(i>0 && Node::Compare_byCountAsc_thenByCodeDesc(*a[i-1], *n) > 0) {
				i--;
			}
			a.insert(a.begin()+i, n);
		}
	}

	// Decompressor can't handle trees of <2 leafs. Make sure there are at least 2.
	if (a.size() < 2)
	{
		if (a.size() == 0) {
			stats[1] = 1;
			stats[2] = 1;
		} 
		else {
			int c = a[0]->Chr;
			if (c-1 >= 0)
				stats[c - 1] = 1;
			else
				stats[c + 1] = 1;
		}
		goto reconstruct_a;
	}

	while (a.size() > 1)
	{
		nodes[nq++] = Node(a[0], a[1], a[0]->Count + a[1]->Count);
		Node* n = &nodes[nq-1];
		a.erase(a.begin()+0, a.begin()+2);
		// высокие ноды кладём в конец, так дерево будет лучше сбалансировано
		size_t i = a.size(); // insertion sort
		while(i>0 && Node::Compare_byCountAsc_thenByCodeDesc(*a[i-1], *n) > 0) {
			i--;
		}
		a.insert(a.begin()+i, n);
	}

	vector<int> code_lengths(stats.size());
	function<void(Node*, int)> traverse = [&](Node* n, int depth) {
		if (n->IsLeaf) {
			code_lengths[n->Chr] = depth;
		} else {
			traverse(n->L, depth + 1);
			traverse(n->R, depth + 1);
		}
	};
	traverse(a[0], 0);

	// Если получились слишком длинные коды, используем спецалгоритм, ограничивающий длину кодов
	int maxlen = 0;
	for(int d : code_lengths) {
		maxlen = max(maxlen, d);
	}
	if (maxlen > codeLenLimit) {
		code_lengths = MakeCodeTable_LimitedLen(stats, codeLenLimit);
	}

	vector<string> ct = rip_tree_builder::MakeCodes(code_lengths);

	for (int c = 0; c < (int)stats.size(); c++) {
		if (ct[c].length() != code_lengths[c])
			throw exception(); // rip_tree_builder выдал другие длины кодов
	}

	return ct;
};

void Write_trees(BitStream& output, vector<string>& charTree, vector<string>& offsetTree)
{
	// Temp table
	vector<int> lengthsTable;
	for (int c = 0; c < 0x120; c++) lengthsTable.push_back((int)charTree[c].length());
	for (int c = 0; c < 32; c++) lengthsTable.push_back((int)offsetTree[c].length());
	if (lengthsTable[EOS_CODE] == 0) {
		throw exception("EOS code has zero length");
	}

	// The use of 0x11 code saves 5.8 bytes on average.
	// Try several freq for code 17 to find optimal one.
	int x = lengthsTable.back();
	lengthsTable.push_back(x); // fictitious item 0x140
	lengthsTable.push_back(x); // fictitious item 0x141
	vector<string> bestTempTree;
	BitStream bestBits;
	for (int freq17 = 0; freq17 <= 0x200; freq17 = max(1, freq17 * 2))
	{
		vector<int> tempstats(18);
		for (int l : lengthsTable) {
			if (l > 15) {
				throw exception("too long codes");
			}
			tempstats[l]++;
		}
		tempstats[16] = 1; // repeat + end of table char
		tempstats[17] = freq17; // repeat
		vector<string> tempTree = MakeCodeTable(tempstats, MAX_CODE_LEN);

		BitStream bits;
		bool eosWritten = false;
		for (int i = 0; i < 0x140; ) {
			if (i > 0 && lengthsTable[i] == lengthsTable[i-1] && lengthsTable[i+1] == lengthsTable[i-1]) {
				if (i + 2 >= 0x140)	{
					bits.Write(tempTree[16]); // repeat + eos
					i += 2;
					eosWritten = true;
				}
				else if (tempTree[17] != "" && tempTree[17].length() < 2 * tempTree[lengthsTable[i]].length()) {
					bits.Write(tempTree[17]); // repeat
					i += 2;
				}
				else {
					bits.Write(tempTree[lengthsTable[i++]]);
				}
			}
			else {
				bits.Write(tempTree[lengthsTable[i++]]);
			}
		}
		if (!eosWritten) {
			bits.Write(tempTree[16]); // repeat + eos
		}

		if (bestBits.Count == 0 || bits.Count < bestBits.Count) {
			bestBits = bits;
			bestTempTree = tempTree;
		}
	}

	for (int c = 0; c < 18; c++) {
		output.Write(4, (int)bestTempTree[c].length());
	}

	for (size_t j = 0; j < bestBits.Count; j++) {
		output.Write(1, bestBits[j]);
	}
}

struct PackResult
{
	vector<byte> PackedData;
	vector<int> ResultCodeStats, ResultOffsetStats;
};

struct Op
{
	ushort Len;
	ushort Offset; // 0 for reuse offset; #FFxx means Char

	int Char() {
		return Offset >= 0xFF00 ? (byte)Offset : -1;
	}

	Op() {
		memset(this, 0, sizeof(*this));
	}

	Op(int len, int offset) {
		if (len < 1 || len > MAX_LEN + MAX_LEN_CORRECTION)
			throw exception();
		if (offset < 0 || offset > MAX_OFFSET)
			throw exception();
		Len = (ushort)len;
		Offset = (ushort)offset;
	}
	
	Op(byte chr) {
		Len = 0;
		Offset = 0xFF00 + chr;
	}
};

PackResult Pack(const vector<byte>& data_x2, int dataSize_x1, vector<string>& cncTree, vector<string>& offsetTree, int pass,
				int longestMatch, const vector<int*>& resPool, const vector<Op*>& resOpPool)
{
	BitStream output;
	Write_trees(output, cncTree, offsetTree);

	// precalc chars encoding cost
	vector<byte> charCost(256);
	for (int c = 0; c < 256; c++) {
		charCost[c] = cncTree[c] != "" ? (byte)cncTree[c].length() : 255;
	}

	// precalc len encoding codes and cost
	vector<string> lenCodes(49152);
	vector<sbyte> lenCost(49152);
	vector<byte> lenChr(49152);
	for (int a = 1; a < 32; a++) {
		int len = a; // 01..1f
		if (a >= 5)	{
			int t = (a + 0xFB + 2) & 0xFF; // 02..1c
			int extrabits = t / 2; // 01..0e
			for (int i = 0; i < (1 << extrabits); i++) {
				len = 2 + (t & 1); // 2..3
				len = (len << extrabits) + i;
				len++; // 5..49152
				lenChr[len - 1] = 255;
				lenCost[len - 1] = -1;
				if (cncTree[256 + a] != "") { // allowed len?
					string code = cncTree[256 + a];
					for (int j = extrabits - 1; j >= 0; j--) {
						code += i >> j & 1;
					}
					lenCodes[len - 1] = code;
					lenCost[len - 1] = (sbyte)code.length();
					lenChr[len - 1] = (byte)a;
				}
			}
		}
		else { // 1..4
			lenChr[len - 1] = 255;
			lenCost[len - 1] = -1;
			if (cncTree[256 + a] != "") {
				lenCodes[len - 1] = cncTree[256 + a];
				lenCost[len - 1] = (sbyte)cncTree[256 + a].length();
				lenChr[len - 1] = (byte)a;
			}
		}
	}

	// Skip list of allowed match lengths
	vector<ushort> nextAllowedLen(49152);
	int k = 0xffff;
	for (int i = 49152; i > 0; i--) {
		nextAllowedLen[i - 1] = (ushort)k;
		if (lenCost[i - 1] >= 0)
			k = i;
	}

	// precalc offset encoding codes and cost
	vector<string> offsetCodes(49152);
	vector<sbyte> offsetCost(49152);
	vector<byte> offsetChr(49152);
	for (int a = 1; a < 32; a++) {
		int offset = a; // 01..1f
		if (a >= 5) {
			int t = (a + 0xFB + 2) & 0xFF; // 02..1c
			int extrabits = t / 2; // 01..0e
			for (int i = 0; i < (1 << extrabits); i++) {
				offset = 2 + (t & 1); // 2..3
				offset = (offset << extrabits) + i;
				offset++; // 5..49152
				offsetChr[offset - 1] = 255;
				offsetCost[offset - 1] = -1;
				if (offsetTree[a] != "") { // allowed offset?
					string code = offsetTree[a];
					for (int j = extrabits - 1; j >= 0; j--) {
						code += i >> j & 1;
					}
					offsetCodes[offset - 1] = code;
					offsetCost[offset - 1] = (sbyte)code.length();
					offsetChr[offset - 1] = (byte)a;
				}
			}
		}
		else { // 1..4
			offsetChr[offset - 1] = 255;
			offsetCost[offset - 1] = -1;
			if (offsetTree[a] != "") {
				offsetCodes[offset - 1] = offsetTree[a];
				offsetCost[offset - 1] = (sbyte)offsetTree[a].length();
				offsetChr[offset - 1] = (byte)a;
			}
		}
	}

	int reuseOffsetCost = offsetTree[0] != "" ? (int)offsetTree[0].length() : -1;


	const byte* data = data_x2.data();
	const int N = dataSize_x1;

	vector<int> matchLen(N);
	if (N * 2 != data_x2.size()) {
		throw exception();
	}
	const byte* zdata = data_x2.data();

	vector<int*> res(N + 1);
	vector<Op*> resOp(N + 1);
	{
		// For position N result is all zeroes
		int max_last_offset = max(0, min(N-1, MAX_OFFSET));
		//res[N] = vector<int>(max_last_offset+1);
		res[N] = resPool[N];
		memset(res[N], 0, sizeof(int) * (max_last_offset+1));
	}

	// Find solution (using dynamic programming and exhaustive search)

	double lastreport = 0;
	for (int de = N - 1; de >= 0; de--)
	{
		double now = clock() / (double)CLOCKS_PER_SEC;
		if (now < lastreport || now > lastreport + 0.03) { // limit to 30 reports/sec
			printf("\rpass %d: %d     ", pass, de);
			lastreport = now;
		}

		ZFunction(zdata, de, N, matchLen.data());
		//for (int hl = de - 1; hl >= 0; hl--)
		//	matchLen[N - (de - hl)] = min(matchLen[N - (de - hl)], N - de);

		// suppose we are in position de=5. Then last_offset <= 4
		int max_last_offset = max(0, min(de - 1, MAX_OFFSET));

		//resOp[de] = vector<Op>(max_last_offset+1);
		resOp[de] = resOpPool[de];
		memset(resOp[de], 0, sizeof(Op) * (max_last_offset + 1));

		//res[de] = vector<int>(max_last_offset+1);
		if (de + longestMatch + 1 <= N) {
			// we don't need results for position de+LM+1 anymore, so let's reuse memory
			res[de] = res[de + longestMatch + 1];
			res[de + longestMatch + 1] = NULL;
		} else {
			res[de] = resPool[de];
		}
		memset(res[de], 0x7f, sizeof(int) * (max_last_offset + 1)); // 0x7f7f7f7f, almost infinity :)

		// Try normal offsets
		int longmatch_best = INT_MAX; // infinity
		Op longmatch_bestOp = Op();
		{
			int lenLowerBound = 0;
			int min_hl = max(0, de - MAX_OFFSET);
			for (int hl = de - 1; hl >= min_hl; hl--)
			{
				int offset = de - hl;
				int ocost = offsetCost[offset - 1];
				if (ocost >= 0) // allowed offset?
				{
					int lenCorrection = offset >= 0x100 ? 1 : 0;
					int max_len = min(min(matchLen[N - offset], N-de), MAX_LEN + lenCorrection);
					int len = 1 + lenCorrection;
					len = max(len, lenLowerBound -10);
					int encodedLen = len - lenCorrection;
					for (; len <= max_len; len = lenCorrection + (encodedLen = nextAllowedLen[encodedLen - 1])) {
						int u = lenCost[encodedLen - 1];
						if (u >= 0) { // allowed len?
							int t = ocost + u + res[de + len][offset];
							if (t < longmatch_best) {
								longmatch_best = t; longmatch_bestOp = Op(len, offset);
							}
						}
					}
					lenLowerBound = max_len;
				}
			}
		}

		// Compute result for different last_offset's
		int lenLowerBound = 0;
		for (int last_offset = 0; last_offset <= max_last_offset; last_offset++)
		{
			// try insert byte
			int best = charCost[data[de]] + res[de + 1][last_offset];
			Op bestOp = Op(data[de]);

			// Try reuse offset
			if (0 != last_offset && reuseOffsetCost >= 0)
			{
				int o = last_offset;
				int minl = max(1, lenLowerBound -10);
				int maxl = min(min(matchLen[N - o], N-de), MAX_LEN);
				for (int l = minl; l <= maxl; l = nextAllowedLen[l - 1]) {
					int u = lenCost[l - 1];
					if (u >= 0) // allowed len?
					{
						int t = reuseOffsetCost + u + res[de + l][o];
						if (t <= best) {
							best = t; bestOp = Op(l, 0);
						}
					}
				}
				lenLowerBound = maxl;
			}

			// Try match with new offset
			if (longmatch_best < best) {
				best = longmatch_best; bestOp = longmatch_bestOp;
			}

			res[de][last_offset] = best;
			resOp[de][last_offset] = bestOp;
		} // last_offset

	} // de

	// Emit compressed data
	vector<int> realCodeStats(0x120);
	vector<int> realOffsetStats(32);
	int lastOffset = 0;
	int resultSizeInBits = res[0][lastOffset];
	size_t writtenBits0 = output.Count;
	{
		int de;
		for (de = 0; de < N; )
		{
			Op op = resOp[de][lastOffset];
			if (op.Char() >= 0) {
				if (op.Char() != data[de]) {
					throw exception();
				}
				output.Write(cncTree[op.Char()]);
				de++;
				realCodeStats[op.Char()]++;
			}
			else if (op.Offset == 0) // use lastOffset
			{
				if (lastOffset <= 0 || lastOffset > de || op.Len < 1) {
					throw exception();
				}
				if (lenCodes[op.Len - 1] == "" || offsetTree[0] == "") {
					throw exception();
				}
				output.Write(lenCodes[op.Len - 1]);
				output.Write(offsetTree[0]);
				de += op.Len;
				realCodeStats[lenChr[op.Len - 1] + 256]++;
				realOffsetStats[0]++;
			}
			else {
				int lenCorrection = op.Offset >= 0x100 ? 1 : 0;
				int encodedLen = op.Len - lenCorrection;
				if (encodedLen < 1) {
					throw exception();
				}
				if (lenCodes[encodedLen - 1] == "" || offsetCodes[op.Offset - 1] == "") {
					throw exception();
				}
				output.Write(lenCodes[encodedLen - 1]);
				output.Write(offsetCodes[op.Offset - 1]);
				de += op.Len;
				realCodeStats[lenChr[encodedLen - 1] + 256]++;
				realOffsetStats[offsetChr[op.Offset - 1]]++;
				lastOffset = op.Offset;
			}
		}
		if (de != dataSize_x1) {
			throw exception();
		}
	}

	size_t writtenBits = output.Count - writtenBits0;
	if (writtenBits != resultSizeInBits) {
		throw exception("writtenBits != resultSizeInBits");
	}

	output.Write(cncTree[EOS_CODE]);
	realCodeStats[EOS_CODE]++;

	//printf("\rpass %d:      ", pass);
	//printf("\rpass %d: ", pass);

	PackResult pr;
	pr.PackedData = output.Data;
	pr.ResultCodeStats = realCodeStats;
	pr.ResultOffsetStats = realOffsetStats;
	return pr;
}

PackResult Pack(const vector<byte>& data_x2, int size, vector<int>& charstats, vector<int>& offsetStats, int pass,
				int longestMatch, vector<int*>& resPool, vector<Op*>& resOpPool)
{
	if (charstats[EOS_CODE] != 1) {
		throw exception("EOS code stats is invalid");
	}
	vector<string> codetable = MakeCodeTable(charstats, MAX_CODE_LEN);
	vector<string> offsettable = MakeCodeTable(offsetStats, MAX_CODE_LEN);
	return Pack(data_x2, size, codetable, offsettable, pass, 
				longestMatch, resPool, resOpPool);
}

vector<byte> Pack(const vector<byte>& data_x2, int dataSize, std::function<void(vector<byte>&)>& saveResult)
{
	if (dataSize * 2 != data_x2.size()) {
		throw exception();
	}
	const byte* data = data_x2.data();

	// Precalc required amount of memory
	int longestMatch = FindLongestMatch(data, dataSize);
	longestMatch = min(longestMatch, MAX_LEN + MAX_LEN_CORRECTION);
	longestMatch = max(1, longestMatch); // at least one position ahead is required for insert_char op
	//printf("Longest match: %d\n", longestMatch);

	int64 memreq = 0;
	for (int pos = dataSize; pos >= 0; pos--) {
		if (pos < dataSize) {
			memreq += (int64)sizeof(Op) * min(pos, MAX_OFFSET + 1);
		}
		if (pos + longestMatch + 1 > dataSize) {
			memreq += (int64)sizeof(int) * min(pos, MAX_OFFSET + 1);
		}
	}
	printf("Allocating %I64d MiB of memory\n", memreq >> 20);

	// Preallocate memory, so we don't waste time if it's not enough
	vector<int*> resPool(dataSize + 1);
	vector<Op*> resOpPool(dataSize + 1);
	for (int pos = dataSize; pos >= 0; pos--)
	{
		int max_lastoffset = max(0, min(pos-1, MAX_OFFSET));
		if (pos < dataSize) {
			if (!(resOpPool[pos] = new Op[max_lastoffset + 1]))
				throw bad_alloc();
		}
		if (pos + longestMatch + 1 <= dataSize) {
			// when processing position X, we don't need results for position X+LM+1 and above,
			// so let's reuse memory by reusing array from position X+LM+1
		} else {
			if (!(resPool[pos] = new int[max_lastoffset + 1]))
				throw bad_alloc();
		}
	}

	auto t0 = clock();

	// Сделаем начальную таблицу кодов и используем её для первого прохода.
	// Цель первого прогона: узнать примерно, какие коды насколько полезны.

	int maxreplen = dataSize - 1;
	int maxoffset = dataSize - 1;

	vector<int> charstats(0x120);
	for (int i=0; i < dataSize; i++) {
		charstats[data[i]] = 1;
	}
	charstats[256] = 1; // EOS_CODE
	for (int lenc = 1; lenc <= 31; lenc++) {
		 charstats[257 + lenc - 1] = (maxreplen >= range_floor[lenc] ? 1 : 0);
	}

	vector<int> offsetStats(32);
	offsetStats[0] = 1; // reuse offset
	for(int c = 1; c <= 31; c++) {
		offsetStats[c] = (maxoffset >= range_floor[c] ? 1 : 0);
	}

	PackResult best;
	int retries = 0;
	for(int pass = 1; ; pass++)
	{
		PackResult packres = Pack(data_x2, dataSize, charstats, offsetStats, pass, longestMatch, resPool, resOpPool);
		printf("\rpass %d: %d bytes\n", pass, packres.PackedData.size());
		if (pass == 1 || packres.PackedData.size() < best.PackedData.size())
		{
			best = packres;
			saveResult(packres.PackedData);
			retries = 0;
		}
		else
		{
			if (++retries >= 2)
				break; // результат не улучшался на протяжении 3-х проходов
		}

		// Возьмём статистику кодов от последнего прохода и используем её для следующего
		charstats = packres.ResultCodeStats;
		offsetStats = packres.ResultOffsetStats;
	}

	timeSpent = double(clock() - t0) / CLOCKS_PER_SEC;

	// Release memory
	for (int pos = dataSize; pos >= 0; pos--) {
		if (resPool[pos]) delete[] resPool[pos];
		if (resOpPool[pos]) delete[] resOpPool[pos];
	}

	return best.PackedData;
}

vector<byte> compress(const vector<byte>& data_x2, int datasize_x1, std::function<void(vector<byte>&)>& saveResult)
{
	timeSpent = -1;
	auto r = Pack(data_x2, datasize_x1, saveResult);
	return r;
};


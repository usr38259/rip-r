/*
 * Copyright © 2022 Eugene Larchenko <el6345@gmail.com>. All rights reserved.
 * See the attached LICENSE file.
 */

#include "longest_match.h"
#include "types.h"
#include <exception>

// Find max L, 0<=L<=n-1, such that there exists i,j such that a{i..i+L-1} == a{j..j+L-1}
int FindLongestMatch(const byte* a, int n)
{
	if (n < 2) {
		return 0;
	}

	// Using binary search to find L.
	// Using sliding hash and a hashtable to find matches of given len.
	// Collisions are possible, in which case the result will be greater than it should, which is acceptable.

	const int M1 = 0x000ffffd; // some small prime
	const int M2 = 0x007ffff1; // some bigger prime
	if (n > M1 / 2) {
		throw std::exception("file is too large"); // need larger hashtable
	}
	int* hashtable = new int[M1];
	if (!hashtable) {
		throw std::bad_alloc();
	}
	int min = 0, max = n-1;
	while (min < max)
	{
		int L = (min + max + 1) / 2;

		bool foundMatch = false;
		for (int j = 0; j < M1; j++) hashtable[j] = -1;
		int i;
		int hash1 = 0, c1 = 1;
		int hash2 = 0, c2 = 1;
		for (i = 0; i < L; i++)
		{
			hash1 = (hash1 * 256 + a[i]) % M1;
			c1 = (c1 * 256) % M1;
			hash2 = (hash2 * 256 + a[i]) % M2;
			c2 = (c2 * 256) % M2;
		}
		hashtable[hash1] = hash2;
		for (; i < n; i++)
		{
			hash1 = ((hash1 * 256 - c1 * a[i-L] + a[i]) % M1 + M1) % M1;
			hash2 = ((hash2 * 256 - c2 * a[i-L] + a[i]) % M2 + M2) % M2;

			int j;
			for (j = hash1; hashtable[j] >= 0; j = (j+1) % M1) {
				if (hashtable[j] == hash2) {
					foundMatch = true;
					goto k0;
				}
			}

			hashtable[j] = hash2;
		}

	k0:
		if (foundMatch) {
			min = L;
		} else {
			max = L-1;
		}
	}

	delete[] hashtable;

	if (min != max)
		throw std::exception(); // impossible

	return min;
}

#pragma once

#include "types.h"
#include <vector>
#include <functional>

using namespace std;

extern double timeSpent;
extern vector<byte> compress(const vector<byte>& data_x2, int datasize_x1, std::function<void(vector<byte>&)>& saveResult);

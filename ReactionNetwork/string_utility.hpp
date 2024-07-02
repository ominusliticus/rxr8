#pragma once

#include "print.hpp"
#include <cstddef>
#include <string>

#define SPACE 32
#define TAB   9

inline bool
is_space(char c)
{
	return (int)c == SPACE || (int)c == TAB;
}

[[nodiscard]] inline std::vector<std::string>
split_string(std::string const& line)
{
	std::vector<std::string> split_string;
	std::size_t              len = line.size();
	std::size_t              pos{ 0 };

	bool in_word{ false };
	for (std::size_t n{ 0 }; n < len; ++n)
	{
		if (!in_word)
		{
			if (!is_space(line[n]))
			{
				if (n == len - 1) { split_string.push_back(line.substr(n)); }
				in_word = true;
				pos     = n;
			}
		}
		else
		{
			if (is_space(line[n]))
			{
				split_string.push_back(line.substr(pos, n - pos));
				in_word = false;
			}
			else if (n == len - 1) split_string.push_back(line.substr(pos));
		}
	}

	return split_string;
}
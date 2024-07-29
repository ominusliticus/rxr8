#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>

template<typename Stream, typename T>
Stream&
operator<<(Stream& stream, std::vector<T> const& vec)
{
	stream << "{ ";
	for (auto const& x : vec)
		stream << x << ", ";
	stream << "}" << std::endl;
	return stream;
}

template<typename Stream, typename Key, typename Value>
Stream&
operator<<(Stream& stream, std::unordered_map<Key, Value> const& umap)
{
	stream << "{ ";
	for (auto const& [key, value] : umap)
		stream << key << ": " << value << ", ";
	stream << "}" << std::endl;
	return stream;
}

template<typename... Args>
void
print(Args&&... args)
{
	((std::cout << std::forward<Args>(args) << " "), ...);
	std::cout << std::endl;
}

template<char delim, typename... Args>
void
print_delim(Args&&... args)
{
	((std::cout << std::forward<Args>(args) << delim), ...);
	std::cout << std::endl;
}

template<typename Stream, typename... Args>
void
fprint(Stream& stream, Args&&... args)
{
	((stream << std::forward<Args>(args) << " "), ...);
	stream << std::endl;
}

template<char delim, typename Stream, typename... Args>
void
fprint_delim(Stream& stream, Args&&... args)
{
	((stream << std::forward<Args>(args) << delim), ...);
	stream << std::endl;
}
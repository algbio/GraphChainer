#ifndef SERIALIZER_HPP
#define SERIALIZER_HPP

#include <string>
#include <vector> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <unordered_map>

class Serializer {
public:
	std::fstream file;
public:
	Serializer() {}
	Serializer(const std::string &filename) : file(filename, std::ios::out | std::ofstream::binary) {}

	void close() {
		file.close();
	}

	template<typename T>
	Serializer& operator &(T &t) {
		return (*this) << t;
	}

	Serializer& operator <<(int t) {
		file.write((char *)&t, sizeof(t));
		return *this;
	}

	Serializer& operator <<(long long t) {
		file.write((char *)&t, sizeof(t));
		return *this;
	}

	Serializer& operator <<(size_t t) {
		file.write((char *)&t, sizeof(t));
		return *this;
	}

	Serializer& operator <<(std::string &s) {
		(*this) << s.length();
		file.write(s.c_str(), s.length());
		return *this;
	}

	template<typename T>
	Serializer& operator <<(std::vector<T> &vec) {
		(*this) << vec.size();
		for (T &t : vec)
			(*this) << t;
		return *this;
	}

	template<typename T1, typename T2>
	Serializer& operator <<(std::pair<T1, T2> &t) {
		return (*this) << t.first << t.second;
	}

	template<typename T1, typename T2>
	Serializer& operator <<(std::unordered_map<T1, T2> &vec) {
		(*this) << vec.size();
		for (std::pair<T1, T2> &t : vec)
			(*this) << t;
		return *this;
	}

	template<typename T>
	Serializer& operator <<(T &t) {
		return t.Serialize(*this);
	}
};

class StringSerializer {
public:
	std::ostringstream ss;
	bool first;
	std::string sep, parL, parR;
public:
	StringSerializer(const std::string &sep = ", ", const std::string &parL = "{ ", const std::string &parR = " }") : first(true), sep(sep), parL(parL), parR(parR) {}
	std::string str() {
		return ss.str();
	}

	template<typename T>
	static std::string to_string(const T &t, const std::string &sep = ", ", const std::string &parL = "{ ", const std::string &parR = " }") {
		StringSerializer s(sep, parL, parR);
		s << t;
		return s.str();
	}

	template<typename T>
	static std::string to_string(T &t, const std::string &sep = ", ", const std::string &parL = "{ ", const std::string &parR = " }") {
		StringSerializer s(sep, parL, parR);
		s << t;
		return s.str();
	}

	void putsep() {
		if (!first) ss << sep;
		first = false;
	}

	template<typename T>
	StringSerializer& operator &(const T &t) {
		return (*this) << t;
	}

	StringSerializer& operator <<(int t) {
		putsep();
		ss << t;
		return *this;
	}

	StringSerializer& operator <<(long long t) {
		putsep();
		ss << t;
		return *this;
	}

	StringSerializer& operator <<(size_t t) {
		putsep();
		ss << t;
		return *this;
	}

	StringSerializer& operator <<(const std::string &s) {
		putsep();
		ss << "\"" << s << "\"";
		return *this;
	}

	template<typename T>
	StringSerializer& operator <<(const std::vector<T> &vec) {
		putsep();
		ss << parL;
		first = true;
		for (const T &t : vec)
			(*this) << t;
		ss << parR;
		return *this;
	}

	template<typename T1, typename T2>
	StringSerializer& operator <<(const std::pair<T1, T2> &t) {
		putsep();
		ss << parL;
		first = true;
		(*this) << t.first;
		(*this) << t.second;
		ss << parR;
		return *this;
	}

	template<typename T1, typename T2>
	StringSerializer& operator <<(const std::unordered_map<T1, T2> &vec) {
		putsep();
		ss << parL;
		first = true;
		for (const std::pair<T1, T2> &t : vec)
			(*this) << t;
		ss << parR;
		return *this;
	}

	template<typename T>
	StringSerializer& operator <<(const T &t) {
		putsep();
		ss << parL;
		first = true;
		auto &ar = t.Serialize(*this);
		ss << parR;
		return ar;
	}
};

class Deserializer {
public:
	std::fstream file;
public:
	Deserializer() {}
	Deserializer(const std::string &filename) : file(filename, std::ios::in | std::ofstream::binary) {}

	void close() {
		file.close();
	}

	template<typename T>
	Deserializer& operator >>(T &t) {
		return t.Serialize(*this);
	}

	template<typename T>
	Deserializer& operator &(T &t) {
		return (*this) >> t;
	}

	Deserializer& operator >>(int &t) {
		file.read((char *)&t, sizeof(t));
		return *this;
	}

	Deserializer& operator >>(long long &t) {
		file.read((char *)&t, sizeof(t));
		return *this;
	}

	Deserializer& operator >>(size_t &t) {
		file.read((char *)&t, sizeof(t));
		return *this;
	}

	Deserializer& operator >>(std::string &s) {
		size_t n;
		(*this) >> n;
		s.resize(n);
		if (n > 0)
			file.read(&s[0], n);
		return *this;
	}

	template<typename T>
	Deserializer& operator >>(std::vector<T> &vec) {
		size_t n;
		(*this) >> n;
		vec.resize(n);
		for (size_t i = 0; i < n; i++)
			(*this) >> vec[i];
		return *this;
	}

	template<typename T1, typename T2>
	Deserializer& operator >>(std::pair<T1, T2> &t) {
		return (*this) >> t.first >> t.second;
	}

	template<typename T1, typename T2>
	Deserializer& operator >>(std::unordered_map<T1, T2> &vec) {
		vec.clear();
		size_t n;
		(*this) >> n;
		vec.reserve(n);
		T1 t1;
		T2 t2;
		for (size_t i = 0; i < n; i++) {
			(*this) >> t1 >> t2;
			vec[t1] = t2;
		}
		return *this;
	}
};

#endif
/// \file util/StringBuilder.h
/// \brief Utility class for quick stream to string conversion.
#ifndef UTIL_STRINGBUILDER_H
#define UTIL_STRINGBUILDER_H

#include <string>
#include <sstream>
#include <string.h>

namespace util
{

/// Utility class for quick stream to string conversion.
class StringBuilder
{
	std::ostringstream m_oss;

public:
	/// Stream operator
	template<typename T> StringBuilder& operator << (const T& t)
	{
		m_oss << t;
		return *this;
	}

	/// String conversion operator.
	operator std::string() const
	{
		return m_oss.str();
	}

	/// c_str() conversion function.
	const char* c_str() const
	{
		return m_oss.str().c_str();
	}
};

/// \brief istringstream wrapper that can scan for literal (const) strings too.

/// Utility class to extend standard istringstream functionality with
/// the ability to match string literals (and set fail bit accordingly).
/// This enables the ability to, like sscanf(), also take literal formatting into account.
/// E.g.
///    StringScanner ss("foo40_bar2");
///    ss >> "foo" >> i >> "bar" >> j;
///    return !ss.fail()
///
class StringScanner
{
	std::istringstream m_iss;
	char m_delim;

public:

	/// Constructors
	StringScanner( std::ios_base::openmode which = std::ios_base::in ) :
	  m_iss(which),
	  m_delim('\0')
	{ }

	StringScanner( const std::string & str, std::ios_base::openmode which = std::ios_base::in ) :
	  m_iss(str, which),
	  m_delim('\0')
	{ }

	/// Get the std::istreamstream object
	inline
	std::istringstream & iss()
	{
		return m_iss;
	}

	/// \brief Set a delimiter character, for scanning std::string objects ('\0' to unset it
	///
	/// Example:
	///   ss("foo hello|world"); ss << s1; ss.delim('|') << s2 << s3;
	/// \see StringScanner& operator >> (std::string &str)
	inline
	StringScanner& delim(char delim) {
		m_delim = delim;
		return *this;
	}

	/// Stream operator
	template<typename T>
	inline
	StringScanner& operator >> (T& t)
	{
		m_iss >> t;
		return *this;
	}

	inline
	StringScanner& operator >> (std::string &str)
	{
		if (m_delim == '\0') {
			// use normal stream operator
			m_iss >> str;
		}
		else {
			// get string upto delimiter is found,
			// discard delimiter
			std::getline(m_iss, str, m_delim);
		}
		return *this;
	}

	/// Special case for string literals (C-style string with \0 delimiter),
	/// constant string should match literally, and be discarded
	/// otherwise fail bit will be set.
	StringScanner& operator >> (const char *token)
	{
		// get token length and create buffer
		size_t l = strlen(token);
		char *buffer = new char[l];

		// test if token is found on stream
		m_iss.read(buffer, l);
		if (strncmp(buffer, token, l) != 0)
		{
			// token could not be found,
			// set fail flag
			m_iss.clear(std::ios::failbit);
		}

		// cleanup buffer
		delete[] buffer;

		return *this;
	}
	
	/// Encapsulation of istringstream::fail()
	inline
	bool fail ( ) const {
		return m_iss.fail();
	}

	/// Encapsulation of istringstream::bad()
	inline
	bool bad ( ) const {
		return m_iss.bad();
	}
	
	/// Encapsulation of istringstream::good()
	inline
	bool good ( ) const {
		return m_iss.good();
	}

	/// Encapsulation of istringstream::eof()
	inline
	bool eof ( ) const {
		return m_iss.eof();
	}

	/// Encapsulation of istringstream::rdstate()
	inline
	std::ios::iostate rdstate ( ) const {
		return m_iss.rdstate();
	}

	/// Encapsulation of istringstream::clear()
	inline
	void clear( std::ios::iostate state = std::ios::goodbit ) {
		m_iss.clear(state);
	}

	/// Encapsulation of istringstream::read()
	inline
	StringScanner& read ( char* s, std::streamsize n )
	{
		m_iss.read(s, n);
		return *this;
	}

	/// Encapsulation of istringstream::peek()
	inline
	int peek ( ) {
		return m_iss.peek();
	}

};

// use this define for convenient short notation
#define UTIL_STR(args) (util::StringBuilder() << args)

/// \brief Convenient shorthand for util::StringScanner usage.
///
/// Example:
///   if (UTIL_SCAN("foo4_bar2", "foo" >> i >> "_bar" >> j))
///   { ... }
///
#define UTIL_SCAN(str, args) (!(util::StringScanner(str) >> args).fail())


} // namespace util

#endif // UTIL_UTIL_STRINGBUILDER_H

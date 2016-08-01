/// \file util/StringFormat.h
/// \brief Utility class for more convenient printf/scanf type functionality.
///
/// Example
///    StringFormatter format("image_{06d:frame}_{d:cam}.png");
///
///    format.set("frame", 12).set("cam", 3);
///    cout << format.print() << endl;
///
///    format.scan("image_001234_6.png");
///    int cam = format.get("cam");
///

#ifndef VBL_UTIL_STRINGFORMAT_H
#define VBL_UTIL_STRINGFORMAT_H

#include <string>
#include <sstream>
#include "StringBuilder.h"

namespace util
{

class ParamBase
{
public:

	std::string m_name;

	ParamBase(const std::string &name) :
		m_name ( name )
	{ }

	ParamBase(const ParamBase & other ) :
		m_name (other.m_name)
	{ }

	virtual ~ParamBase() {}

	virtual ParamBase * copy() const = 0;
	virtual void print(util::StringBuilder &builder) = 0;
	virtual void scan(util::StringScanner &scanner) = 0;
};

template<typename T>
class Param : public ParamBase
{
public:
	T m_value;

	Param<T>(const std::string &name) :
		ParamBase(name)
	{
	}

	Param<T>(const Param<T> & other ) :
		ParamBase(other.m_name), m_value(other.m_value)
	{ }

	ParamBase * copy() const
	{
		return new Param<T>(*this);
	}

	void print(util::StringBuilder &builder)
	{
		builder << m_value;
	}

	void scan(util::StringScanner &scanner)
	{
		scanner >> m_value;
	}
};

template<>
class Param<std::string> : public ParamBase
{
public:
	std::string m_value;

	Param<std::string>(const std::string &name) : 
		ParamBase(name)
	{ }

	ParamBase * copy() const
	{
		return new Param<std::string>(*this);
	}

	void print(util::StringBuilder &builder)
	{
		builder << m_value;
	}

	void scan(util::StringScanner &scanner)
	{
		if (m_name.empty())
		{
			// this is a constant string, match the scanned string exactly
			// (StringScanner knows by casting to C-string)
			scanner >> m_value.c_str();
		}
		else
		{
			scanner >> m_value;
		}
	}
};

template<>
class Param<int> : public ParamBase
{
public:
	int m_value;
	int m_width;
	char m_fill;

	Param<int>(const std::string &name, int width = 0, char fill = '0') : 
		ParamBase(name), m_width (width), m_fill(fill)
	{ }

	ParamBase * copy() const
	{
		return new Param<int>(*this);
	}

	void print(util::StringBuilder &builder)
	{
		builder << std::setw(m_width) << std::setfill(m_fill) << m_value;
	}

	void scan(util::StringScanner &scanner)
	{
		scanner >> m_value;
	}
};

template<>
class Param<float> : public ParamBase
{
public:
	float m_value;
	int m_precision;

	Param<float>(const std::string &name, int precision = 5) : 
		ParamBase(name), m_precision (precision)
	{ }

	ParamBase * copy() const
	{
		return new Param<float>(*this);
	}

	void print(util::StringBuilder &builder)
	{
		builder << std::setprecision(m_precision) << m_value;
	}

	void scan(util::StringScanner &scanner)
	{
		scanner >> m_value;
	}
};

/// \brief Utility class for more convenient printf/scanf type functionality.

/// Example
///    StringFormatter format("image_{06d:frame}_{d:cam}.png");
///
///    format.set("frame", 12).set("cam", 3);
///    cout << format.print() << endl;
///
///    format.scan("image_001234_6.png");
///    int cam = format.get("cam");
///
class StringFormatter
{
	std::vector<ParamBase *> m_params;
	std::map<std::string, ParamBase *> m_namemap;

public:
	/// Default constructor
	StringFormatter()
	{ }

	/// Copy constructor
	StringFormatter(const StringFormatter & other)
	{
		std::map<std::string, ParamBase *>::const_iterator end = other.m_namemap.end();

		for(size_t j = 0; j < other.m_params.size(); ++j)
		{
			ParamBase * pBase = other.m_params[j]->copy(); 
			m_params.push_back(pBase);

			if (other.m_namemap.find(pBase->m_name) != end)
			{
				m_namemap[pBase->m_name] = pBase;
			}
		}	
	}

	/// Destructor, deletes parameter objects
	~StringFormatter()
	{
		for(size_t j = 0; j < m_params.size(); ++j)
		{
			delete m_params[j];
		}
	}

	/// Initialize with a format string
	StringFormatter(const std::string &format)
	{
		parseFormat(format);
	}

	/// Stream operator
	template<typename T> StringFormatter& operator << (const Param<T> & param)
	{
		Param<T> * pParam = new Param<T>(param);
		m_params.push_back(pParam);
		m_namemap[param.m_name] = pParam;

		return *this;
	}

	StringFormatter& operator << (const char * fixed_value)
	{
		return fixed(fixed_value);
	}

	StringFormatter& fixed(std::string fixed_value)
	{
		Param<std::string> * pParam = new Param<std::string>("");
		pParam->m_value = fixed_value;
		m_params.push_back(pParam);

		// do not store this in m_namemap, thus cannot be changed
		return *this;
	}

	// set value
	template<typename T>
	StringFormatter & set(const std::string name, T value)
	{
		ParamBase * pBase = m_namemap[name];
		if (pBase)
		{
			// update value
			Param<T> * pParam = dynamic_cast<Param<T> *>(pBase);
			
			VBL_ASSERT(pParam != NULL, "Named parameter not of given type");

			pParam->m_value = value;
		}
		{
			// no such named parameter
		}

		return *this;
	}

	// get value
	template<typename T>
	StringFormatter & get(const std::string name, T& value)
	{
		ParamBase * pBase = m_namemap[name];
		if (pBase)
		{
			Param<T> * pParam = dynamic_cast<Param<T> *>(pBase);
			value = pParam->m_value;
		}
		{
			// no such named parameter
		}

		return *this;
	}

	template<typename T>
	T get(const std::string name)
	{
		T value;
		get(name, value);
		return value;
	}

	// sprintf / sscanf style functions

	/// printf-like construction of a string
	std::string print() const
	{
		util::StringBuilder builder;

		for(size_t j = 0; j < m_params.size(); ++j)
		{
			m_params[j]->print(builder);
		}

		return builder;
	}

	/// scanf-like parsing of a string
	bool scan(const std::string &input, bool testEOF = true)
	{
		util::StringScanner scanner(input);

		for(size_t j = 0; j < m_params.size(); ++j)
		{
			m_params[j]->scan(scanner);
		}

		if (scanner.fail())
			return false;

		if (!testEOF)
			return true; // no need to test if full input was consumed

		// test that we fully finished scanning (reached end of stream)
		scanner.peek(); // forces EOF bit to be set if we really reached the end
		return scanner.eof();
	}

	/// String conversion operator.
	operator std::string() const
	{
		return print();
	}

	/// c_str() conversion function.
	const char* c_str() const
	{
		return print().c_str();
	}

	/// Parse format string
	void parseFormat(const std::string & format)
	{
		StringScanner scanner(format);

		const char DELIM_OPEN = '{';
		const char DELIM_CLOSE = '}';
		std::string token;

		scanner.delim(DELIM_OPEN);
		while(!scanner.eof())
		{
			scanner >> token;
			fixed(token); // fixed string matching

			if (scanner.eof())
			{
				// last fixed-string token found
				break;
			}

			// find "{arg}" argument by looking for closing tag
			scanner.delim(DELIM_CLOSE);
			scanner >> token;
			parseFormatArgument(token);

			// reset
			scanner.delim(DELIM_OPEN);
		}
		fixed("\0");
	}

	void parseFormatArgument(const std::string & format)
	{
		std::string type, name;
		char ctype;
		
		int width = 0, precision = 5;
		char fill = ' ';

		StringScanner scanner(format);
		scanner.delim(':');
		scanner >> type >> name;

		ctype = type.at(type.length()-1); // last character

		// create part
		switch (ctype)
		{
			case 'd':

				if (type.at(0) == '0') { fill = '0'; } // check for "06" type format
				VBL_SCAN(type, width); // parse starting number as width

				*this << Param<int>(name, width, fill);
				break;

			case 'f':

				VBL_SCAN(type, precision); // parse starting number as precision

				*this << Param<float>(name, precision);
				break;

			case 's':

				*this << Param<std::string>(name);
				break;

			default:
				throw std::exception("Incorrect format type specified");
				
		}
	}

};

// use this define for convenient short notation
#define VBL_STRING_FORMAT(args) (util::StringFormatter() << args)
#define VBL_FORMAT(args) (util::StringFormatter(args))

} // namespace util

#endif // VBL_UTIL_STRINGFORMAT_H

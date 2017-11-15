#ifndef EXCEPTIONDEF
#define EXCEPTIONDEF


//  *****************
//  *  Exception class  *
//  *****************


//  Exceptions thrown by the class of vectors are dealt with by this class

//  "throw" throws two strings
//  The first string summarises the error
//  The second string describes the error
//  A constructor is written to construct objects from these two strings

//  the function DebugPrint() prints details of the error




#include <string>

class Exception
{
public:
  std::string problem, summary;
  Exception(std::string sum, std::string problem);
  void DebugPrint();
};

#endif

#ifndef APPENDERFACTORY_H_
#define APPENDERFACTORY_H_

namespace proteomics_sdk
{
class CAppender;
class CAppenderFactory
{
public:
	CAppenderFactory();
	virtual ~CAppenderFactory();
	CAppender *GetAppender(LogAppenderType eType)const;
	
};

}

#endif /*APPENDERFACTORY_H_*/

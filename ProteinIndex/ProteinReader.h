#ifndef PROTEINREADER_H_
#define PROTEINREADER_H_
#include <string>
using namespace std;


namespace ProteinIndex
{
class CProteinReader {
public:
	virtual ~CProteinReader() {
	}
	virtual void Open(const std::string &strWorkDir, const std::string & strDBName) = 0;
	virtual void Close() = 0;
	virtual bool GetNext(string& strProtein) = 0;
	virtual size_t GetCurrentID() const = 0;
protected:
	size_t m_tProteinID;
};


class CRandomReader : public CProteinReader {
public:
	virtual ~CRandomReader() {
	}
//	virtual void Open(const string &strFilePath) = 0;
	virtual void Open(const std::string &strWorkDir, const std::string & strDBName) = 0;
	virtual void Close() = 0;

	virtual bool GetNext(string& strProtein) = 0;
	virtual size_t GetCurrentID() const = 0;
	virtual void SetCurrentID(size_t tCurrentID) = 0;
	virtual size_t GetByID( string& , size_t proID ) = 0;
	virtual size_t GetProNum(void) const = 0;
//	virtual long CalcPosInPointerBolck(size_t proID) const = 0;//liyou

//	virtual long GetPos(size_t proID) const = 0;
//	virtual long GetPos(unsigned char uFileID) const = 0;
	
	virtual void ReadPepSQ(string& strDesSQ, size_t tProID, size_t tStartPos, unsigned char cLen)  = 0;
	//virtual void ReadPepSQ(string& strDesSQ, unsigned char cFileID, size_t tStartPos, unsigned char cLen)  = 0;
	virtual void ReadPepSQ_EX(string& strDesSQ, size_t tProID, size_t tStartPos, unsigned char cLen)  = 0;
	//virtual void ReadPepSQ_EX(string& strDesSQ, unsigned char cFileID,size_t tPos, unsigned char cLen) =0;
	
protected:
	virtual inline void SetWorkDir(const char * szDestDir)
	{
		string strDestDir( szDestDir );
		string::size_type idx = 0;
		idx = strDestDir.find('\\', idx);
		while (idx != string::npos ) {
			strDestDir[idx] = '/';
			idx = strDestDir.find('\\', idx);
		}

		idx = strDestDir.find_last_not_of(" ");
		if (strDestDir[idx] != '/' ) 
		{
			m_strWorkDir = strDestDir.substr(0,idx+1) + "/";
		}
		else
			m_strWorkDir = strDestDir;
			return;
	}
	
	virtual inline string GetWorkDir(void) const 
	{
		return m_strWorkDir;
	}
	
	virtual inline void SetDBName(const char * szDBName)
	{
		m_strDBName = szDBName;
	}
	
	virtual inline string GetDBName(void) const
	{
		return m_strDBName;
	}



	string	m_strWorkDir;
	// The name (or path) for protein sequence database.
	string	m_strDBName;
	size_t m_tProNum;
};
}

#endif /*PROTEINREADER_H_*/

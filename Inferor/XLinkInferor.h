#ifndef XLINK_INFEROR_H_
#define XLINK_INFEROR_H_

class CXLinkInferor : public CInferor
{
public:
	CXLinkInferor();
	virtual ~CXLinkInferor();
	
	virtual bool Init(string strOption, time_t tmStartTime = 0);
	virtual void Run(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier);
	
private:
	void _WritePfd2PlabelConfigFile(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier);
	void _WritePfd2PXBuildConfigFile(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier);
	void _WritePfd2PBuildConfigFile(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier);
	//members:
	CCondition m_Condition;

	string m_strOption;
	string m_strWorkDir;
	
	string m_strSpectraPath;
	
	time_t m_tmStartTime;

	CTrace *m_pTrace;
};

#endif /*XLINK_INFEROR_H_*/

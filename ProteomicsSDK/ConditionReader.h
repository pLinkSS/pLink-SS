#ifndef CONDITIONREADER_H_
#define CONDITIONREADER_H_

namespace proteomics_sdk
{
class CCondition;

/* CConditionReader: read parameter files */
class CConditionReader {
public:
	CConditionReader(const string &strOption, const string &strWorkDir);
	virtual ~CConditionReader();
	
	void Read();
	
	static void Test(const CCondition & c);
protected:
	void _LoadDatabase(void);
	void _LoadEnzyme(void);
	void _LoadModify(void);
	void _LoadXLinker(void);
	void _LoadFlow(void);

	void _LoadSimpleScore(void); /* add for xlink open flow */
	void _LoadIons(void);
	void _LoadSpectra(void);
	void _LoadCluster(void);
	void _LoadMSearchInfo(void);
	void _LoadMetaDat(); // add by czhou
public:
	string m_strOption; /* absolute path of parameter file *.pfind */
	string m_strWorkDir; /* installation path of the software */
	
	CCondition m_Condition; /* variable to store searching parameters */
}; /* end of class CConditionReader */
} /* end of namespace proteomics_sdk */

#endif /* end of CONDITIONREADER_H_*/

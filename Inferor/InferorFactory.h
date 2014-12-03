#ifndef INFERORFACTORY_H_
#define INFERORFACTORY_H_

class CInferorFactory
{
public:
	CInferorFactory();
	virtual ~CInferorFactory();
	CInferor * GetInferor(INFEROR_TYPE eType);
};

#endif /*INFERORFACTORY_H_*/

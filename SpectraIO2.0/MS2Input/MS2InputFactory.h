/*
 * MS2InputFactory.h
 *
 *  Created on: 2009-7-7
 *      Author: fan
 */

#ifndef MS2INPUTFACTORY_H_
#define MS2INPUTFACTORY_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"



class CMS2InputFactory {
public:
	CMS2InputFactory();
	virtual ~CMS2InputFactory();

	CMS2Input * GetImporter(MS2FormatType eType);
	
};


#endif /* MS2INPUTFACTORY_H_ */

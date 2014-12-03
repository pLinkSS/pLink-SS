/*
 * MS2InputFactory.h
 *
 *  Created on: 2009-7-7
 *      Author: fan
 */

#ifndef MS2OUTPUTFACTORY_H_
#define MS2OUTPUTFACTORY_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"



class CMS2OutputFactory {
public:
	CMS2OutputFactory();
	virtual ~CMS2OutputFactory();

	CMS2Output * GetImporter(MS2FormatType eType);
};


#endif /* MS2INPUTFACTORY_H_ */

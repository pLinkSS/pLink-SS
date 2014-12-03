#include "../include/sdk.h"
#include "../include/interface.h"
#include "ACInfer.h"
#include "ProteinInferFactory.h"

using namespace proteomics_sdk;
#define _DEBUG_WPWANG_MPI

#ifdef _DEBUG_WPWANG_MPI
#endif  

CProteinInfer * CProteinInferFactory::GetProteinInfer(ProteinInferType eProteinInferType)
{
	switch(eProteinInferType)
	{
	case PROTEIN_INFER_DEFAULT:
		return new CACInfer();
	default:
		return new CACInfer();
	}
}

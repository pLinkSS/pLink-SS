#pragma once


using namespace proteomics_sdk;



namespace proteomics_sdk
{
class CProteinInferFactory
{
public:
	CProteinInfer * GetProteinInfer(ProteinInferType eProteinInferType);
};
}


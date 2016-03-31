#pragma once

typedef struct tagGrabCutParameters 
{
	double			dLambda;
	unsigned int	nNumFGMMs;
	unsigned int	nNumBGMMs;
	unsigned int	nNoiseRadius;
}GrabCutParameters;

typedef struct tagTextureParameters 
{
	unsigned int	nTensorScale;
	double			dTensorEsp;
	double			dTensorTVPower;
	unsigned int	nTensorNosteps;
	double			dTensorSigma;
	double			dTensorStepsize;
	unsigned int	nLowDim;
}TextureParameters;

class Options
{
public:
	Options(void);
	~Options(void);

	GrabCutParameters m_grabcut;
	TextureParameters m_texture;
};

inline Options::Options()
{
	m_grabcut.dLambda = 50;
	m_grabcut.nNumBGMMs = 5;
	m_grabcut.nNumFGMMs = 5;
	m_grabcut.nNoiseRadius = 1;
	m_texture.nTensorScale = 2;
	m_texture.dTensorSigma = 0;
	m_texture.dTensorEsp = 0.001;
	m_texture.dTensorTVPower = 0.6;
	m_texture.dTensorStepsize = 5000;
	m_texture.nTensorNosteps = 2;
	m_texture.nLowDim = 2;
}

inline Options::~Options()
{
}
#include "mmFilter.h"
#include <mmreg.h>
#include <ks.h>
#include <ksmedia.h>

#include "mex.h"

// unfortunately to detect memory leaks, we have to put code in each .cpp file
#if defined(_DEBUG) && defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif

mmFilter::mmFilter() :
	CSource(NAME("mmFilter"), NULL, CLSID_mmFilter)
{ }

HRESULT mmFilter::addVideo(char** frames, double* times, int nr, int width, int height)
{
	HRESULT hr;
	videoSources.add(new videoSource(this, frames, times, nr, width, height,&hr));
	return hr;
}

HRESULT mmFilter::addAudio(char** frames, double* times, int nr, int* lens, int nrChannels, int rate)
{
	HRESULT hr;
	audioSources.add(new audioSource(this, frames, times, nr, lens, nrChannels, rate, &hr));
	return hr;
}

mmFilter::~mmFilter()
{
	_RPT1(_CRT_WARN,"Freeing video Sources %d\n",videoSources.size());
	for (int i=0; i<videoSources.size(); i++) delete videoSources.at(i);
	_RPT1(_CRT_WARN,"Freeing audio Sources %d\n",audioSources.size());
	for (int i=0; i<audioSources.size(); i++) delete audioSources.at(i);
}


mmSource::mmSource(CSource *pParent, char** frames, double* times, int	nr, HRESULT* phr) :
	CSourceStream(NAME("mmSourceStream"), phr, pParent, L"mmSourcePin"),
	CSourceSeeking(NAME("mmSourceSeeking"), (IPin*) this, phr, &m_cSharedState),
	currentFrame(0),
	done(0),
	stopGraph(1)
{
	this->frames = frames;
	this->times = times;
	this->nr = nr;
	m_dRateSeeking = 1;
	m_rtStart = 0;

	sleepEvent = CreateEvent(NULL, TRUE, TRUE, "sleepEvent"); 
}

mmSource::~mmSource()
{ }

STDMETHODIMP mmSource::NonDelegatingQueryInterface( REFIID riid, void ** ppv )
{
	if(	riid ==	IID_IMediaSeeking ) 
	{
	return CSourceSeeking::NonDelegatingQueryInterface( riid, ppv );
	}
	return CSourceStream::NonDelegatingQueryInterface(riid, ppv);
}

void mmSource::UpdateFromSeek()
{
	if (ThreadExists())	
	{
		DeliverBeginFlush();
		Stop();
		DeliverEndFlush();
		Run();
	}
}

HRESULT	mmSource::OnThreadStartPlay()
{
	_RPT3(_CRT_WARN,"OnThreadStartPlay %d %d %f\n",(int)m_rtStart, (int)m_rtStop, (float)m_dRateSeeking);
	return DeliverNewSegment(m_rtStart, m_rtStop, m_dRateSeeking);
}

HRESULT	mmSource::ChangeStart( )
{
	_RPT0(_CRT_WARN,"ChangeStart\n");
	{
		CAutoLock lock(CSourceSeeking::m_pLock);
		currentFrame = 0;
	}
	UpdateFromSeek();
	return S_OK;
}

HRESULT	mmSource::ChangeStop( )
{
	_RPT0(_CRT_WARN,"ChangeStop\n");
	{
		CAutoLock lock(CSourceSeeking::m_pLock);
		if (currentFrame < nr && times[currentFrame]*10000000 < m_rtStop) return S_OK;
	}
	UpdateFromSeek();
	return S_OK;
}

HRESULT	mmSource::OnThreadCreate()
{
	_RPT0(_CRT_WARN,"OnThreadCreate\n");

	return NOERROR;
}

void mmSource::SetStopGraph(int stopGraph)
{
//	CAutoLock cAutoLockShared(&m_cSharedState);
	this->stopGraph = stopGraph;
}



audioSource::audioSource(CSource* filt, char**	frames,	double* times, int nr, int* lens, int nrChannels, int rate, HRESULT* phr) :
	mmSource(filt,frames,times,nr,phr)
{
	this->lens = lens;
	this->nrChannels = nrChannels;
	this->rate = rate;
	
	ResetFrames(frames,times,nr,lens);
}

void audioSource::ResetFrames(char** frames, double* times, int nr, int* lens)
{
	_RPT0(_CRT_WARN,"ResetFrames\n");
	CAutoLock cAutoLock(m_pFilter->pStateLock());

	this->frames = frames;
	this->times = times;
	this->nr = nr;
	this->lens = lens;
	currentFrame = 0;
	done = 0;
	
	maxlen = 0;
	for (int i=0;i<nr;i++) if (maxlen < lens[i]) maxlen = lens[i];

	m_rtStop = times[nr-1]*10000000+((LONGLONG)lens[nr-1])*10000000/rate;
}

HRESULT	audioSource::FillBuffer(IMediaSample *pMediaSample)
{	
	_RPT2(_CRT_WARN,"FillBuffer %d %d\n",currentFrame,nr);

	long lDataLen = pMediaSample->GetSize();;
	{
		CAutoLock cAutoLockShared(&m_cSharedState);
		
		if (currentFrame >= nr || times[currentFrame]*10000000 > m_rtStop)
		{
			_RPT0(_CRT_WARN,"a stopping\n");
			done=1;
			
			if (stopGraph) return S_FALSE;
			else {
/*
				pMediaSample->SetActualDataLength(0);
				REFERENCE_TIME rtStart,	rtStop;
				
				rtStart	= times[currentFrame-1]*10000000;
				rtStop = m_rtStop;
				pMediaSample->SetTime(&rtStart,	&rtStop);
			
				_RPT0(_CRT_WARN,"a Sleeping \n");
//				while (currentFrame >= nr || times[currentFrame]*10000000 > m_rtStop) Sleep(1000);
				Sleep(100);
				return NOERROR;
*/
				ResetEvent(sleepEvent);
				while (WAIT_OBJECT_0 != WaitForSingleObject(sleepEvent,INFINITE));
			}
		}
		
		double* dData = (double*)frames[currentFrame];

		if (subtype == MEDIASUBTYPE_PCM)
		{
			short *pData;
			pMediaSample->GetPointer((BYTE**)&pData);
			
			for (int i=lens[currentFrame]*nrChannels-1;i>=0;i--) pData[i] = min((1<<15)-1,dData[i]*(1<<15));
		} else { // FLOAT format
			float *pData;
			pMediaSample->GetPointer((BYTE**)&pData);

			for (int i=lens[currentFrame]*nrChannels-1;i>=0;i--) pData[i] = dData[i];
		}

		REFERENCE_TIME rtStart,	rtStop;
		
		rtStart	= times[currentFrame]*10000000;
		rtStop	= rtStart + lens[currentFrame]/(rate*10000000.0);

		pMediaSample->SetActualDataLength(lens[currentFrame]*wordSize);
		_RPT4(_CRT_WARN,"a SetTime %d %d   %d %d\n",(int)(rtStart>>32),(int)rtStart,(int)(rtStop>>32),(int)rtStop);
		pMediaSample->SetTime(&rtStart,	&rtStop);
	
		currentFrame++;
	}

	pMediaSample->SetSyncPoint(TRUE);

	return NOERROR;
}

HRESULT	audioSource::GetMediaType(int iPosition, CMediaType *pmt)
{
	_RPT0(_CRT_WARN,"GetMediaType\n");
	CAutoLock cAutoLock(m_pFilter->pStateLock());
	_RPT1(_CRT_WARN,"iPosition %d\n",iPosition);
	if (iPosition < 0) return E_INVALIDARG;
//	if (iPosition > 1) return VFW_S_NO_MORE_ITEMS;
	if (iPosition > 0) return VFW_S_NO_MORE_ITEMS;

	if (!pmt) return E_INVALIDARG;

	WAVEFORMATEX *fmt;
	if (nrChannels > 2 || rate > 48000)
	{
		WAVEFORMATEXTENSIBLE* fmte = (WAVEFORMATEXTENSIBLE*) pmt->AllocFormatBuffer(sizeof(WAVEFORMATEXTENSIBLE));
		if (NULL == fmte) return E_OUTOFMEMORY;
		ZeroMemory(fmte, sizeof(WAVEFORMATEXTENSIBLE));
		fmte->Format.cbSize = sizeof(WAVEFORMATEXTENSIBLE)-sizeof(WAVEFORMATEX);
		fmte->Samples.wValidBitsPerSample = iPosition==1?32:16;
		fmte->SubFormat = iPosition==1?KSDATAFORMAT_SUBTYPE_IEEE_FLOAT:KSDATAFORMAT_SUBTYPE_PCM;
		
//		if (nrChannels == 6) fmte->dwChannelMask = KSAUDIO_SPEAKER_5POINT1; // assume 6 channels means 5.1
		
		fmt = (WAVEFORMATEX*) fmte;
		fmt->wFormatTag	= WAVE_FORMAT_EXTENSIBLE;
	} else {
		fmt = (WAVEFORMATEX*) pmt->AllocFormatBuffer(sizeof(WAVEFORMATEX));
		if (NULL == fmt) return E_OUTOFMEMORY;
		ZeroMemory(fmt, sizeof(WAVEFORMATEX));
		fmt->cbSize = 0;
		fmt->wFormatTag	= iPosition==1?WAVE_FORMAT_IEEE_FLOAT:WAVE_FORMAT_PCM;
	}

	fmt->nChannels = nrChannels; 
	fmt->nSamplesPerSec = rate;
	fmt->wBitsPerSample = iPosition==1?32:16;
	fmt->nBlockAlign = fmt->nChannels*fmt->wBitsPerSample/8;
	fmt->nAvgBytesPerSec = fmt->nSamplesPerSec*fmt->nBlockAlign;

	pmt->SetType(&MEDIATYPE_Audio);
	pmt->SetFormatType(&FORMAT_WaveFormatEx);
	pmt->SetSubtype(iPosition==1?&MEDIASUBTYPE_IEEE_FLOAT:&MEDIASUBTYPE_PCM);
	pmt->SetTemporalCompression(FALSE);
	pmt->SetSampleSize(fmt->nBlockAlign); //maxlen*

	return NOERROR;
}

HRESULT	audioSource::CheckMediaType(const CMediaType *pMediaType)
{
	_RPT0(_CRT_WARN,"CheckMediaType\n");
	CAutoLock cAutoLock(m_pFilter->pStateLock());

	if (!pMediaType) return E_INVALIDARG;

	if (*pMediaType->Type() != MEDIATYPE_Audio) return E_INVALIDARG;
	if (*pMediaType->Subtype() != MEDIASUBTYPE_PCM && *pMediaType->Subtype() != MEDIASUBTYPE_IEEE_FLOAT) return E_INVALIDARG;

	WAVEFORMATEX *fmt = (WAVEFORMATEX*) pMediaType->Format();

	if (fmt == NULL) return E_INVALIDARG;

	if ((fmt->nSamplesPerSec != rate) || (fmt->nChannels != nrChannels)) return E_INVALIDARG;

	return S_OK;
}

HRESULT	audioSource::DecideBufferSize(IMemAllocator *pAlloc,ALLOCATOR_PROPERTIES *pProperties)
{
	_RPT0(_CRT_WARN,"DecideBufferSize\n");
	CAutoLock cAutoLock(m_pFilter->pStateLock());

	if (!pAlloc || !pProperties) return E_INVALIDARG;

	ALLOCATOR_PROPERTIES Request, Actual;

	if (pProperties->cbAlign <= 0) Request.cbAlign = 1;
	else Request.cbAlign = pProperties->cbAlign;
	
	if (pProperties->cbBuffer < maxlen*wordSize) Request.cbBuffer = maxlen*wordSize;
	else Request.cbBuffer = pProperties->cbBuffer;

	if (pProperties->cbPrefix < 0) Request.cbPrefix = 0;
	else Request.cbPrefix = pProperties->cbPrefix;

	if (pProperties->cBuffers < 2) Request.cBuffers = 2;
	else Request.cBuffers = pProperties->cBuffers;

	HRESULT hr = pAlloc->SetProperties(&Request,&Actual);
	if (SUCCEEDED(hr)) 
	{
		hr = pAlloc->Commit();
	}

	return hr;
}

HRESULT	audioSource::SetMediaType(const CMediaType *pMediaType)
{
	_RPT0(_CRT_WARN,"SetMediaType\n");
	CAutoLock cAutoLock(m_pFilter->pStateLock());

	if (!pMediaType) return E_INVALIDARG;

	subtype=*pMediaType->Subtype();

	_RPT1(_CRT_WARN,"SetMediaType %s\n",(subtype == MEDIASUBTYPE_PCM)?"PCM":"FLOAT");
	
	WAVEFORMATEX *fmt = (WAVEFORMATEX*) pMediaType->Format();
	wordSize = fmt->nBlockAlign;
	
	HRESULT hr = CSourceStream::SetMediaType(pMediaType);

	_RPT2(_CRT_WARN,"SetMediaType hr %d %d\n",hr,FAILED(hr));
	
	return hr;
}





videoSource::videoSource(CSource* filt, char**	frames,	double* times, int nr, int width, int height, HRESULT* phr) :
	mmSource(filt,frames,times,nr,phr)
{
	this->width = width;
	this->height = height;
	
	ResetFrames(frames,times,nr);
}

void videoSource::ResetFrames(char** frames, double* times, int nr)
{
	_RPT0(_CRT_WARN,"ResetFrames\n");
	CAutoLock cAutoLock(m_pFilter->pStateLock());

	this->frames = frames;
	this->times = times;
	this->nr = nr;
	currentFrame = 0;
	
	m_rtStop = times[nr-1]*10000000;
}

HRESULT	videoSource::FillBuffer(IMediaSample *pMediaSample)
{
//mexWarnMsgTxt("FillBuffer");
	_RPT2(_CRT_WARN,"FillBuffer %d %d\n",currentFrame,nr);

	BYTE *pData;
	long lDataLen = pMediaSample->GetSize();;
	pMediaSample->GetPointer(&pData);
	
	if (lDataLen < height*scanwidth*3) return E_INVALIDARG;
	
	{
		CAutoLock cAutoLockShared(&m_cSharedState);
		if (currentFrame >= nr || times[currentFrame]*10000000 > m_rtStop) 
		{
			_RPT0(_CRT_WARN,"v stopping\n");
			done=1;
			
			if (stopGraph) return S_FALSE;
			else {


/*				pMediaSample->SetActualDataLength(0);
				REFERENCE_TIME rtStart,	rtStop;
				
				rtStart	= times[currentFrame-1]*10000000;
				rtStop = m_rtStop;
				pMediaSample->SetTime(&rtStart,	&rtStop);
			
				_RPT0(_CRT_WARN,"v Sleeping \n");
//				while (currentFrame >= nr || times[currentFrame]*10000000 > m_rtStop) Sleep(30);
				Sleep(100);
				return NOERROR;
*/
				ResetEvent(sleepEvent);
				while (WAIT_OBJECT_0 != WaitForSingleObject(sleepEvent,INFINITE));
			}
		}
	
		//height,width,color => flip color,width, flip height
		BYTE* frame = (BYTE*)frames[currentFrame];

		for (int c=0; c<3; c++)
		{
			int c1 = height*width*c;
			int c2 = 2-c;
			for (int w=0; w<width; w++)
			{
				int w1 = height*w;
				int w2 = 3*w;
				for (int h=0; h<height; h++)
				{
					pData[c2+w2+3*scanwidth*(height-1-h)] = frame[c1+w1+h];
				}
			}
		}
	
		REFERENCE_TIME rtStart,	rtStop;
		
		rtStart	= times[currentFrame]*10000000;
		rtStop	= (currentFrame<nr-1)?times[currentFrame+1]*10000000-1:rtStart;
	
		pMediaSample->SetActualDataLength(height*width*3);
		_RPT4(_CRT_WARN,"v SetTime %d %d   %d %d\n",(int)(rtStart>>32),(int)rtStart,(int)(rtStop>>32),(int)rtStop);
		pMediaSample->SetTime(&rtStart,	&rtStop);
	
		currentFrame++;
	}

	pMediaSample->SetSyncPoint(TRUE);
	
	return NOERROR;
}

HRESULT	videoSource::GetMediaType(int iPosition, CMediaType *pmt)
{
	CAutoLock cAutoLock(m_pFilter->pStateLock());

	if (iPosition < 0) return E_INVALIDARG;
	if (iPosition > 0) return VFW_S_NO_MORE_ITEMS;

	VIDEOINFO *pvi = (VIDEOINFO*) pmt->AllocFormatBuffer(sizeof(VIDEOINFO));
	if (NULL == pvi) return(E_OUTOFMEMORY);

	ZeroMemory(pvi, sizeof(VIDEOINFO));
	pvi->bmiHeader.biCompression = BI_RGB;
	pvi->bmiHeader.biBitCount	 = 24;
	pvi->bmiHeader.biSize	= sizeof(BITMAPINFOHEADER);
	pvi->bmiHeader.biWidth	= width;
	pvi->bmiHeader.biHeight	= height;
	pvi->bmiHeader.biPlanes	= 1;
	pvi->bmiHeader.biSizeImage = GetBitmapSize(&pvi->bmiHeader);
	pvi->bmiHeader.biClrImportant   = 0;
	
	scanwidth = pvi->bmiHeader.biSizeImage/3/height;

	SetRectEmpty(&(pvi->rcSource)); // we want the whole image area rendered.
	SetRectEmpty(&(pvi->rcTarget)); // no particular destination rectangle

	double totalTime = (times[nr-1]-times[0])/(nr-1);

	pvi->AvgTimePerFrame = 10000000*totalTime;
	_RPT1(_CRT_WARN,"pvi->AvgTimePerFrame %d\n",pvi->AvgTimePerFrame);
char errstr[100];
sprintf(errstr,"pvi->AvgTimePerFrame %d",pvi->AvgTimePerFrame);

//mexWarnMsgTxt(errstr);

	pmt->SetType(&MEDIATYPE_Video);
	pmt->SetFormatType(&FORMAT_VideoInfo);
	pmt->SetTemporalCompression(TRUE);

	pmt->SetSubtype(&MEDIASUBTYPE_RGB24);
	pmt->SetSampleSize(pvi->bmiHeader.biSizeImage);

	return NOERROR;
}


HRESULT	videoSource::CheckMediaType(const CMediaType *pMediaType)
{
//mexWarnMsgTxt("CheckMediaType");

	CAutoLock cAutoLock(m_pFilter->pStateLock());

	if (!pMediaType) return E_INVALIDARG;

	if (*pMediaType->Type() != MEDIATYPE_Video || !pMediaType->IsFixedSize()) return E_INVALIDARG;
	if (*pMediaType->Subtype() != MEDIASUBTYPE_RGB24) return E_INVALIDARG;

	VIDEOINFO *pvi = (VIDEOINFO*) pMediaType->Format();

	if (pvi == NULL) return E_INVALIDARG;

	if (pvi->bmiHeader.biWidth != width || pvi->bmiHeader.biHeight != height) return E_INVALIDARG;

	return S_OK;
}

HRESULT	videoSource::DecideBufferSize(IMemAllocator* pAlloc,ALLOCATOR_PROPERTIES* pProperties)
{
//mexWarnMsgTxt("DecideBufferSize");

	CAutoLock cAutoLock(m_pFilter->pStateLock());

	VIDEOINFO *pvi = (VIDEOINFO*) m_mt.Format();

	if (!pAlloc || !pProperties || !pvi) return E_INVALIDARG;

	ALLOCATOR_PROPERTIES Request, Actual;

	if (pProperties->cbAlign <= 0) Request.cbAlign = 1;
	else Request.cbAlign = pProperties->cbAlign;
	
	if (pProperties->cbBuffer < pvi->bmiHeader.biSizeImage) Request.cbBuffer = pvi->bmiHeader.biSizeImage;
	else Request.cbBuffer = pProperties->cbBuffer;

	if (pProperties->cbPrefix < 0) Request.cbPrefix = 0;
	else Request.cbPrefix = pProperties->cbPrefix;

	if (pProperties->cBuffers < 2) Request.cBuffers = 2;
	else Request.cBuffers = pProperties->cBuffers;

	HRESULT hr = pAlloc->SetProperties(&Request,&Actual);
	if (SUCCEEDED(hr)) 
	{
		hr = pAlloc->Commit();
//mexWarnMsgTxt("DecideBufferSize Yay");
	}

//mexWarnMsgTxt("DecideBufferSize end");

	return hr;
}

HRESULT	videoSource::SetMediaType(const CMediaType* pMediaType)
{
//mexWarnMsgTxt("SetMediaType");

	CAutoLock cAutoLock(m_pFilter->pStateLock());

	return CSourceStream::SetMediaType(pMediaType);
}



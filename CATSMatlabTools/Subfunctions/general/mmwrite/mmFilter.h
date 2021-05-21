#include <olectl.h>
#include "streams.h"
#include <initguid.h>
#include "mmutils.h"

// { fd501041-8ebe-11ce-8183-00aa00577da2 }
DEFINE_GUID(CLSID_mmFilter,
0xfd501041, 0x8ebe, 0x11ce, 0x81, 0x83, 0x00, 0xaa, 0x00, 0x57, 0x7d, 0xa2);

class mmSource : public CSourceStream, public CSourceSeeking
{
public:
	mmSource(CSource* pParent, char** frames, double* times, int nr, HRESULT* phr);
	virtual ~mmSource();
	
	STDMETHODIMP NonDelegatingQueryInterface( REFIID riid, void ** ppv );
	
	virtual HRESULT FillBuffer(IMediaSample *pms) = NULL;
	virtual HRESULT DecideBufferSize(IMemAllocator *pIMemAlloc, ALLOCATOR_PROPERTIES *pProperties) = NULL;
	virtual HRESULT SetMediaType(const CMediaType *pMediaType) = NULL;
	virtual HRESULT CheckMediaType(const CMediaType *pMediaType) = NULL;
	virtual HRESULT GetMediaType(int iPosition, CMediaType *pmt) = NULL;
	
	STDMETHODIMP Notify(IBaseFilter * pSender, Quality q) { return NOERROR;};
	
	void SetStopGraph(int stopGraph);

	HRESULT OnThreadCreate(void);
	
private:
	void UpdateFromSeek();
	
protected:
	HRESULT ChangeStart();
	HRESULT ChangeStop();
	HRESULT ChangeRate() { return S_OK; }
	HRESULT OnThreadStartPlay();

	CCritSec m_cSharedState;
	IReferenceClock* pClock;
	REFERENCE_TIME startTime;
	CAMEvent syncEvent;
	DWORD syncToken;

	char** frames;
	double* times;
	int nr;
	int currentFrame;
	int stopGraph;

public:
	int done;

	HANDLE sleepEvent;
};

class audioSource : public mmSource
{
public:
	audioSource(CSource* filt, char** frames, double* times, int nr, int* lens, int nrChannels, int rate, HRESULT* phr);
	
	HRESULT FillBuffer(IMediaSample *pms);
	HRESULT DecideBufferSize(IMemAllocator *pIMemAlloc, ALLOCATOR_PROPERTIES *pProperties);
	HRESULT SetMediaType(const CMediaType *pMediaType);
	HRESULT CheckMediaType(const CMediaType *pMediaType);
	HRESULT GetMediaType(int iPosition, CMediaType *pmt);
	
	void ResetFrames(char** frames, double* times, int nr, int* lens);
private:
	int nrChannels;
	int rate;
	int* lens;
	int maxlen;
	int wordSize;
	GUID subtype;
};

class videoSource : public mmSource
{
public:
	videoSource(CSource* filt, char** frames, double* times, int nr, int width, int height, HRESULT* phr);
	
	HRESULT FillBuffer(IMediaSample *pms);
	HRESULT DecideBufferSize(IMemAllocator *pIMemAlloc, ALLOCATOR_PROPERTIES *pProperties);
	HRESULT SetMediaType(const CMediaType *pMediaType);
	HRESULT CheckMediaType(const CMediaType *pMediaType);
	HRESULT GetMediaType(int iPosition, CMediaType *pmt);

	void ResetFrames(char** frames, double* times, int nr);
private:
	int width, height, scanwidth;
};

class mmFilter : public CSource
{
public:
	mmFilter();
	~mmFilter();

	HRESULT addVideo(char** frames, double* times, int nr, int width, int height);
	HRESULT addAudio(char** frames, double* times, int nr, int* lens, int nrChannels, int rate);

	vector<audioSource*> audioSources;
	vector<videoSource*> videoSources;
};

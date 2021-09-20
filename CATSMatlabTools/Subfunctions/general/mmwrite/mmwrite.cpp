#include "mmFilter.h"
#include <dshowasf.h>
//#include <Dmodshow.h>
//#include <DMOReg.h>
//#include <Wmcodecdsp.h>
#include <wchar.h>

// unfortunately to detect memory leaks, we have to put code in each .cpp file
#if defined(_DEBUG) && defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif
#include <stdlib.h>
#include <crtdbg.h>

#include "mex.h"
#include "mmutils.cpp"

#pragma warning(disable : 4995) // get rid of warning about sprintf and things being deprecated

HRESULT insertEncoder(IGraphBuilder* pGraphBuilder, IBaseFilter* encoder, IBaseFilter* pMux, IPin* outPin)
{
	HRESULT hr;
	
	if (encoder)
	{
		pGraphBuilder->AddFilter(encoder, L"encoder");
		
		IPin* encoderInPin = getInputPin(encoder);
		IPin* encoderOutPin = getOutputPin(encoder);

		if (encoderInPin == NULL || encoderOutPin == NULL) mexErrMsgTxt("Can't get Encoder Pins");
		
		_RPT0(_CRT_WARN,"Connecting Encoder Input Pin\n");
		if (FAILED(hr = pGraphBuilder->Connect(outPin,encoderInPin))) return hr;
		encoderInPin->Release();
		outPin->Release();
		
		outPin = encoderOutPin;
	}

	IPin* Muxpin = NULL;
	IEnumPins* pinList;
	int connected = 0;

	_RPT0(_CRT_WARN,"looping through Mux pins\n");
	//get the input
	if (SUCCEEDED(hr = pMux->EnumPins(&pinList)))
	{
		pinList->Reset();
		while (pinList->Next(1, &Muxpin, NULL) == S_OK && !connected)
		{
			IPin* tmp = 0;
			_RPT0(_CRT_WARN,"testing muxpin dir and connection\n");
			if (getPinInfo(Muxpin).dir == PINDIR_INPUT && Muxpin->ConnectedTo(&tmp)==VFW_E_NOT_CONNECTED)
			{
				_RPT0(_CRT_WARN,"Connecting MUX Input Pin\n");
				if (SUCCEEDED(hr = pGraphBuilder->Connect(outPin,Muxpin)))
				{
					connected = 1;
					_RPT0(_CRT_WARN,"CONNECTED!!!\n");
				}/* else {
					IEnumMediaTypes* pEnum;
					outPin->EnumMediaTypes(&pEnum);
					AM_MEDIA_TYPE* pMediaTypes;
					
					pEnum->Reset();
					while (pEnum->Next(1, &pMediaTypes, NULL) == S_OK && !connected)
					{
						_RPT0(_CRT_WARN,"ConnectDirect\n");
						if (SUCCEEDED(hr = pGraphBuilder->ConnectDirect(outPin,Muxpin,pMediaTypes)))
						{
							connected = 1;
							_RPT0(_CRT_WARN,"CONNECTED!!!\n");							
						}
						MyFreeMediaType(*pMediaTypes);
					}
					pEnum->Release();
				}*/
			} else if (tmp) {
				tmp->Release();
				tmp = NULL;
			}
			Muxpin->Release();
		}
		pinList->Release();
	} else return hr;
	
	_RPT1(_CRT_WARN,"connected: %d\n",connected);
	if (connected) outPin->Release();
	
	return (connected)?NOERROR:VFW_E_NOT_CONNECTED;
}

HRESULT ListEncoders(REFCLSID category, char*** list, int* listlen, char* findname, IBaseFilter** foundFilter)
{
	HRESULT hr;
	ICreateDevEnum *pSysDevEnum = NULL;
	IEnumMoniker *pEnum = NULL;
	IMoniker *pMoniker = NULL;
	IPropertyBag *pPropBag = NULL;
	VARIANT var;
	if (!findname)
	{
		*list = NULL;
		*listlen = 0;
	} else {
		*foundFilter = NULL;
	}

	_RPT0(_CRT_WARN,"ListEncoders... \n");
	
	if (FAILED(hr = CoCreateInstance(CLSID_SystemDeviceEnum, NULL, CLSCTX_INPROC_SERVER, IID_ICreateDevEnum, (void**)&pSysDevEnum))) return hr;
	
	hr = pSysDevEnum->CreateClassEnumerator(category, &pEnum, 0); //CLSID_VideoCompressorCategory
	
	if (hr == S_OK)  // S_FALSE means nothing in this category.
	{
		_RPT0(_CRT_WARN,"while (S_OK == pEnum->Next(1, &pMoniker, NULL))\n");
		while (!(foundFilter && *foundFilter) && S_OK == pEnum->Next(1, &pMoniker, NULL))
		{
			pMoniker->BindToStorage(0, 0, IID_IPropertyBag, (void **)&pPropBag);
			VariantInit(&var);
			hr = pPropBag->Read(L"FriendlyName", &var, 0);
			if (SUCCEEDED(hr))
			{
				char name[1000];
				WideCharToMultiByte(CP_ACP,0,var.bstrVal,-1,name,1000-1,NULL,NULL);

				if (findname) _RPT2(_CRT_WARN,"findname '%s' name '%s'\n",findname,name);
				else _RPT1(_CRT_WARN,"name '%s'\n",name);

				if (findname && strcmp(findname,name)==0)
				{
					if (FAILED(hr = pMoniker->BindToObject(NULL, NULL, IID_IBaseFilter, (void**)foundFilter))) return hr;
					_RPT1(_CRT_WARN,"pMoniker->BindToObject hr %d\n",hr);
				}
				
				if (!findname)
				{
					_RPT0(_CRT_WARN,"list1\n");
					*list = (char**)realloc(*list,(++(*listlen))*sizeof(char*));
					_RPT0(_CRT_WARN,"list2\n");
					char* tmp = new char[strlen(name)+1];
					memcpy(tmp,name,strlen(name)+1);
					(*list)[(*listlen)-1] = tmp;
				}
			}

			VariantClear(&var); 
			pPropBag->Release();
			pMoniker->Release();
		}
	}

	pSysDevEnum->Release();
	pEnum->Release();
	
	return S_OK;
}

mmFilter* mmfilt;
IGraphBuilder* pGraphBuilder;
IMediaEvent* pMediaEventEx;
IMediaControl* pMediaControl;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	HRESULT hr;
	WCHAR uniFilename[MAX_PATH] = L"out.avi";
	char filename[MAX_PATH] = "out.avi";
	IBaseFilter* AudioEncoder = NULL;
	IBaseFilter* VideoEncoder = NULL;
	IBaseFilter *pWriter = NULL;
	IBaseFilter *pAviMux = NULL;
	IBaseFilter *pMux = NULL; // don't release this one.
	int nrVideo = 0;
	int nrAudio = 0;
	int audioQuality = 98;
	int videoQuality = 95;
	int outputHeight = 0, outputWidth = 0;
	int inputHeight = 0, inputWidth = 0;
	char prxFile[MAX_PATH] = "";
	int nrAudioChannels = 0;
	int outputFrameRate = -1;
	double aveFrameRate = -1;
	int nrOutput = 0;
	int buildGraph = 1, stopGraph = 1;

	if (nrhs < 1) mexErrMsgTxt("mmwrite must have at least one parameter");

	for (int i=0; i< nrhs; i++) 
	{
		if (mxIsChar(prhs[i])) {
			char param[MAX_PATH];
			mxGetString(prhs[i],param,MAX_PATH);
			if (strcmp(param,"Initialized")==0) {
				buildGraph = 0;
				if (!pGraphBuilder || !mmfilt)
				{
					mexWarnMsgTxt("mmwrite has not been initialized, ignoring 'Initialized' command.");
					buildGraph = 1;
				}
			} else if (strcmp(param,"Continue")==0) {
				stopGraph = 0;
			}
		}
	}
	
	if (buildGraph)
	{
		// Create the graph builder
		if (FAILED(hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, IID_IGraphBuilder, (void**)&pGraphBuilder))) mexErrMsgTxt(message(hr));

		// make sure everything is back to "normal"; cleanUp() also releases the current pGraphBuilder, so we need to make a new one
		cleanUp(&pGraphBuilder);
		if (FAILED(hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, IID_IGraphBuilder, (void**)&pGraphBuilder))) mexErrMsgTxt(message(hr));

		mmfilt = new mmFilter();
		pGraphBuilder->AddFilter(mmfilt,L"mmfilt");
	}

	for (int i=0; i< nrhs; i++) 
	{
		if (mxIsStruct(prhs[i]))
		{
			if (mxGetFieldNumber(prhs[i],"frames") != -1 && mxGetFieldNumber(prhs[i],"width") != -1 && mxGetFieldNumber(prhs[i],"height") != -1)
			{ // video...
				int nr = mxGetNumberOfElements(prhs[i]);
				for (int n=0;n<nr;n++)
				{
					mxArray* framesStruct = mxGetField(prhs[i],n,"frames");
					
					if (!mxIsStruct(framesStruct)) mexErrMsgTxt("video struct, field frames is not a struct as it should be.");
					
					int nrFrames = mxGetNumberOfElements(framesStruct);

					if (nrFrames==0) continue;
					if (mxGetFieldNumber(prhs[i],"times") == -1) mexErrMsgTxt("video struct, field 'times' is missing.");

					char** frames = (char**)newadd(new char*[nrFrames]);
					double* times = mxGetPr(mxGetField(prhs[i],n,"times"));
					int width = mxGetScalar(mxGetField(prhs[i],n,"width"));
					int height = mxGetScalar(mxGetField(prhs[i],n,"height"));

					aveFrameRate = (nrFrames==1)?((times[0]==0)?30:times[0]):(nrFrames-1)/(times[nrFrames-1]-times[0]);
					
					if (mxGetNumberOfElements(mxGetField(prhs[i],n,"times")) != nrFrames) mexErrMsgTxt("the 'times' vector doesn't have the same number of entries as the frames struct.");
					
					if ((width>>1)<<1 != width) mexWarnMsgTxt("the width of the movie isn't even, this may cause problems.");
					if ((height>>1)<<1 != height) mexWarnMsgTxt("the height of the movie isn't even, this may cause problems.");
					
					for (int f=0;f<nrFrames;f++)
					{
						mxArray* cdata = mxGetField(framesStruct,f,"cdata");
						frames[f] = (char*)mxGetPr(cdata);
						if (mxGetNumberOfElements(cdata) != width*height*3) mexErrMsgTxt("error one of the video frames, isn't of the proper size.");
					}
					
					if (buildGraph)
						mmfilt->addVideo(frames, times, nrFrames, width, height);
					else
						mmfilt->videoSources.at(nrVideo)->ResetFrames(frames, times, nrFrames);
						
					nrVideo++;
					
					inputHeight = height;
					inputWidth = width;
				}
			} else if (mxGetFieldNumber(prhs[i],"data") != -1 && mxGetFieldNumber(prhs[i],"rate") != -1) {
				// audio...
				int nr = mxGetNumberOfElements(prhs[i]);
				for (int n=0;n<nr;n++)
				{
					int rate = mxGetScalar(mxGetField(prhs[i],n,"rate"));
					mxArray* data = mxGetField(prhs[i],n,"data");
					int nrChannels = mxGetN(data);
					int nrSamples = mxGetM(data);
					int nrFrames = nrSamples/rate+1; if ((nrFrames-1)*rate == nrSamples) nrFrames--;

					if (nrFrames==0) continue;
					if (mxGetFieldNumber(prhs[i],"times") == -1) mexErrMsgTxt("audio struct, field 'times' is missing.");

					double** frames = (double**)newadd(new double*[nrFrames]);
					int* lens = (int*)newadd(new int[nrFrames]);
					double* datap = mxGetPr(data);
					double* times = (double*)newadd(new double[nrFrames]);
					double* ftimes = mxGetPr(mxGetField(prhs[i],n,"times"));
					
					if (mxGetNumberOfElements(mxGetField(prhs[i],n,"times")) < 1) mexErrMsgTxt("The 'times' vector must at least specify the start time.");
					
					for (int f=0;f<nrFrames;f++)
					{
						lens[f] = (f==nrFrames-1)?nrSamples-rate*f:rate;
						
						frames[f] = (double*)newadd(new double[nrChannels*lens[f]]);
						for (int j=0;j<nrChannels;j++)
							for (int k=0;k<lens[f];k++)
								frames[f][j+k*nrChannels] = datap[f*rate+k+j*nrSamples];
						times[f] = f+ftimes[0];
					}
					
					if (buildGraph)
						mmfilt->addAudio((char**)frames,times,nrFrames,lens,nrChannels,rate);
					else
						mmfilt->audioSources.at(nrAudio)->ResetFrames((char**)frames,times,nrFrames,lens);
					nrAudio++;
					
					if (nrChannels > nrAudioChannels) nrAudioChannels = nrChannels;
				}
			} else if (mxGetFieldNumber(prhs[i],"videoQuality") != -1 || mxGetFieldNumber(prhs[i],"audioQuality") != -1 || 
					mxGetFieldNumber(prhs[i],"outputHeight") != -1 || mxGetFieldNumber(prhs[i],"outputWidth") != -1 || 
					mxGetFieldNumber(prhs[i],"outputFrameRate") != -1 || mxGetFieldNumber(prhs[i],"prxFile") != -1) {
				int num;
				if ((num = mxGetFieldNumber(prhs[i],"videoQuality")) != -1) videoQuality = mxGetScalar(mxGetFieldByNumber(prhs[i],0,num));
				if ((num = mxGetFieldNumber(prhs[i],"audioQuality")) != -1) audioQuality = mxGetScalar(mxGetFieldByNumber(prhs[i],0,num));
				if ((num = mxGetFieldNumber(prhs[i],"outputHeight")) != -1) outputHeight = mxGetScalar(mxGetFieldByNumber(prhs[i],0,num));
				if ((num = mxGetFieldNumber(prhs[i],"outputWidth")) != -1) outputWidth = mxGetScalar(mxGetFieldByNumber(prhs[i],0,num));
				if ((num = mxGetFieldNumber(prhs[i],"outputFrameRate")) != -1) outputFrameRate = mxGetScalar(mxGetFieldByNumber(prhs[i],0,num));
				if ((num = mxGetFieldNumber(prhs[i],"prxFile")) != -1) mxGetString(mxGetFieldByNumber(prhs[i],0,num),prxFile,MAX_PATH);
			} else if (mxGetFieldNumber(prhs[i],"videoCompressor") != -1 || mxGetFieldNumber(prhs[i],"audioCompressor") != -1) {
				int num;
				char name[1000];
				char errstr[1000] = "Can't find encoder.";
				if ((num = mxGetFieldNumber(prhs[i],"videoCompressor")) != -1)
				{
					mxGetString(mxGetFieldByNumber(prhs[i],0,num),name,1000);
					//sprintf(errstr,"Can not find encoder '%s'.",name);

					_RPT1(_CRT_WARN,"Finding: %s\n",name);
					if (FAILED(hr = ListEncoders(CLSID_VideoCompressorCategory, NULL, NULL, name, &VideoEncoder))) mexErrMsgTxt(message(hr));
					if (!VideoEncoder) mexErrMsgTxt(errstr);
				}
				if ((num = mxGetFieldNumber(prhs[i],"audioCompressor")) != -1)
				{
					mxGetString(mxGetFieldByNumber(prhs[i],0,num),name,1000);
//					sprintf(errstr,"Can not find encoder '%s'.",name);

					_RPT1(_CRT_WARN,"Finding: %s\n",name);
					if (FAILED(hr = ListEncoders(CLSID_AudioCompressorCategory, NULL, NULL, name, &AudioEncoder))) mexErrMsgTxt(message(hr));
					if (!AudioEncoder) mexErrMsgTxt(errstr);
				}
			} else mexErrMsgTxt("mmwrite: unrecognized struct format, it must be a video or audio struct, or a configuration struct.  Type 'help mmwrite' for more details.");
		} else if (mxIsChar(prhs[i])) {
			char param[MAX_PATH];
			mxGetString(prhs[i],param,MAX_PATH);

			if (i==0)
			{
				// only the first parameter can be the filename
				strcpy(filename,param);
				MultiByteToWideChar(CP_ACP, 0, filename, -1, uniFilename, MAX_PATH);
			} else {
				if (strcmp(param,"ListAviVideoEncoders")==0)
				{
					char** list;
					int listlen;
					if (FAILED(hr = ListEncoders(CLSID_VideoCompressorCategory, &list, &listlen, NULL, NULL))) mexErrMsgTxt(message(hr));
					
					if (nlhs > nrOutput)
					{
						plhs[nrOutput] = mxCreateCellMatrix(listlen,1);
						for (int j=0;j<listlen;j++)
						{
							mxSetCell(plhs[nrOutput],j,mxCreateString(list[j]));
							delete [] (list[j]);
						}
						nrOutput++;
					}
					if (list) free(list);
					return;
				} else if (strcmp(param,"ListAviAudioEncoders")==0) {
					char** list;
					int listlen;
					if (FAILED(hr = ListEncoders(CLSID_AudioCompressorCategory, &list, &listlen, NULL, NULL))) mexErrMsgTxt(message(hr));
					
					if (nlhs > nrOutput) 
					{
						plhs[nrOutput] = mxCreateCellMatrix(listlen,1);
						for (int j=0;j<listlen;j++)
						{
							mxSetCell(plhs[nrOutput],j,mxCreateString(list[j]));
							delete [] (list[j]);
						}
						nrOutput++;
					}
					if (list) free(list);
					return;
				} else if (strcmp(param,"Initialized")==0) {
					// ignore...
				} else if (strcmp(param,"Continue")==0) {
					// ignore...
				} else mexErrMsgTxt("mmwrite: Unknown command, options are: 'ListAviVideoEncoders', 'ListAviAudioEncoders', 'Initialized', 'Continue'");
			}
		} else mexErrMsgTxt("mmwrite accepts only structs and strings");
	}

	for (int i=0;i<mmfilt->audioSources.size();i++) {mmfilt->audioSources.at(i)->SetStopGraph(stopGraph); mmfilt->audioSources.at(i)->done=0;}
	for (int i=0;i<mmfilt->videoSources.size();i++) {mmfilt->videoSources.at(i)->SetStopGraph(stopGraph); mmfilt->videoSources.at(i)->done=0;}

	if (buildGraph)
	{
		if (strcmp(filename+(strlen(filename)-3),"avi")==0)
		{
			_RPT0(_CRT_WARN,"AVI\n");
			CoCreateInstance(CLSID_AviDest, 0, CLSCTX_INPROC_SERVER, IID_IBaseFilter, (void **)&pAviMux);
			pGraphBuilder->AddFilter(pAviMux, L"AVI MUX");
			pMux = pAviMux;
			
			CoCreateInstance(CLSID_FileWriter, 0, CLSCTX_INPROC_SERVER, IID_IBaseFilter, (void **)&pWriter);
			pGraphBuilder->AddFilter(pWriter, L"Writer");
			
			IPin* in = getInputPin(pWriter);
			IPin* out = getOutputPin(pMux);
			
			if (FAILED(hr = pGraphBuilder->Connect(out,in))) mexErrMsgTxt(message(hr));
			
			in->Release();
			out->Release();
		} else {
			_RPT1(_CRT_WARN,"ASF/WMV '%s'\n",filename+(strlen(filename)-3));
			CoCreateInstance(CLSID_WMAsfWriter, 0, CLSCTX_INPROC_SERVER, IID_IBaseFilter, (void **)&pWriter);
			pMux = pWriter;
			pGraphBuilder->AddFilter(pWriter, L"Writer");
			
			if (VideoEncoder || AudioEncoder) mexErrMsgTxt("Error ASF/WMV/WMA writter is not compatible with specifying an encoder");
/*
			//Try to instantiate the WMV9 Encoder class so that we can tweek the parameters; however, it can fail and keep running
			int noWMV9 = FALSE;
			noWMV9 = FAILED(CoCreateInstance(CLSID_DMOWrapperFilter, 0, CLSCTX_INPROC_SERVER, IID_IBaseFilter, (void **)&VideoEncoder));

			IDMOWrapperFilter* pWrapperFilter = NULL;

			if (!noWMV9) noWMV9 = FAILED(VideoEncoder->QueryInterface(IID_IDMOWrapperFilter, (void **)&pWrapperFilter));

			if (!noWMV9) noWMV9 = FAILED(pWrapperFilter->Init( CLSID_CWMV9EncMediaObject,DMOCATEGORY_VIDEO_ENCODER ));

			IPropertyBag* pProperty = NULL;

			if (!noWMV9) noWMV9 = FAILED(pWrapperFilter->QueryInterface( IID_IPropertyBag, (void**)&pProperty));

			
			if (!noWMV9) {
				VARIANT Var;

				VariantInit(&Var);
				Var.vt = VT_I4;
				Var.lVal = 1;

				pProperty->Write( g_wszWMVCFullFrameRate, &Var);
			}

			if (noWMV9) VideoEncoder = NULL;
*/

			IConfigAsfWriter* pConfigAsfWriter;
			pMux->QueryInterface(IID_IConfigAsfWriter, (void**)&pConfigAsfWriter);
			_RPT1(_CRT_WARN,"pConfigAsfWriter %d\n",pConfigAsfWriter);
			IWMProfile* pProfile, * pProfile2;
			pConfigAsfWriter->GetCurrentProfile(&pProfile);

			IWMProfileManager* pWMProfileManager = NULL;
			WMCreateProfileManager(&pWMProfileManager);

			DWORD dwLength = 0;
			DWORD dwBytesRead = 0;
			HANDLE hFile = INVALID_HANDLE_VALUE;
			WCHAR* data;

			if (strlen(prxFile) > 0)
			{	
				hFile = CreateFile(prxFile,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
				if(INVALID_HANDLE_VALUE == hFile) mexWarnMsgTxt(message(HRESULT_FROM_WIN32(GetLastError())));
				
				dwLength = GetFileSize(hFile, NULL);
				data = (WCHAR*)new BYTE[dwLength + sizeof(WCHAR)];
				ZeroMemory(data, dwLength + sizeof(WCHAR));
				// there tends to be two bytes at the beginning that will corrupt our stream
				if(!ReadFile(hFile, data, 2, &dwBytesRead, NULL)) mexWarnMsgTxt(message(HRESULT_FROM_WIN32(GetLastError())));
				if(!ReadFile(hFile, data, dwLength-2, &dwBytesRead, NULL)) mexWarnMsgTxt(message(HRESULT_FROM_WIN32(GetLastError())));
				CloseHandle(hFile);
			} else {
				data = new WCHAR[10000];
				WCHAR tmp[1000];
				ZeroMemory(data, 10000*sizeof(WCHAR));
	
				if (outputHeight == 0 || outputWidth == 0)
				{
					outputHeight = inputHeight;
					outputWidth = inputWidth;
				}

				int frameTime = (outputFrameRate == -1)? 10000000/aveFrameRate: 10000000/outputFrameRate;
				
				appendStrW(data,L"<profile version=\"589824\" storageformat=\"1\" name=\"\" description=\"\">");
				appendStrW(data,L"<streamconfig majortype=\"{73647561-0000-0010-8000-00AA00389B71}\" streamnumber=\"1\" streamname=\"Audio Stream\" inputname=\"Audio409\" bitrate=\"128016\" bufferwindow=\"-1\" reliabletransport=\"0\" decodercomplexity=\"\" rfc1766langid=\"en-us\" vbrenabled=\"1\" bitratemax=\"0\" bufferwindowmax=\"0\" vbrquality=\"");
				appendStrW(data,Int2Wstr(audioQuality,tmp,10));
				appendStrW(data,L"\">");
				//VBR Audio 9 16bit
				appendStrW(data,L"<wmmediatype subtype=\"{00000161-0000-0010-8000-00AA00389B71}\" bfixedsizesamples=\"1\" btemporalcompression=\"0\" lsamplesize=\"11889\">");
				appendStrW(data,L"<waveformatex wFormatTag=\"353\" nChannels=\"2\" nSamplesPerSec=\"44100\" nAvgBytesPerSec=\"2147483490\" nBlockAlign=\"11889\" wBitsPerSample=\"16\" codecdata=\"008800000F0000000000\"/>");
				appendStrW(data,L"</wmmediatype>");
				appendStrW(data,L"</streamconfig>");
	
				//VBR Video 9
				appendStrW(data,L"<streamconfig majortype=\"{73646976-0000-0010-8000-00AA00389B71}\" streamnumber=\"2\" streamname=\"Video Stream\" inputname=\"Video409\" bitrate=\"500000\" bufferwindow=\"-1\" reliabletransport=\"0\" decodercomplexity=\"AU\" rfc1766langid=\"en-us\" vbrenabled=\"1\" bitratemax=\"0\" bufferwindowmax=\"0\" vbrquality=\"");
				appendStrW(data,Int2Wstr(videoQuality,tmp,10));
				appendStrW(data,L"\">");
				appendStrW(data,L"<videomediaprops maxkeyframespacing=\"20000000\" quality=\"0\"/>");
				appendStrW(data,L"<wmmediatype subtype=\"{33564D57-0000-0010-8000-00AA00389B71}\" bfixedsizesamples=\"0\" btemporalcompression=\"1\" lsamplesize=\"0\">");
				appendStrW(data,L"<videoinfoheader dwbitrate=\"500000\" dwbiterrorrate=\"0\" avgtimeperframe=\"");
				appendStrW(data,Int2Wstr(frameTime,tmp,10));
				appendStrW(data,L"\">");
				appendStrW(data,L"<rcsource left=\"0\" top=\"0\" right=\"");
				appendStrW(data,Int2Wstr(outputWidth,tmp,10));
				appendStrW(data,L"\" bottom=\"");
				appendStrW(data,Int2Wstr(outputHeight,tmp,10));
				appendStrW(data,L"\"/>");
				appendStrW(data,L"<rctarget left=\"0\" top=\"0\" right=\"");
				appendStrW(data,Int2Wstr(outputWidth,tmp,10));
				appendStrW(data,L"\" bottom=\"");
				appendStrW(data,Int2Wstr(outputHeight,tmp,10));
				appendStrW(data,L"\"/>");
				appendStrW(data,L"<bitmapinfoheader biplanes=\"1\" bibitcount=\"24\" bicompression=\"WMV3\" bisizeimage=\"0\" bixpelspermeter=\"0\" biypelspermeter=\"0\" biclrused=\"0\" biclrimportant=\"0\" biwidth=\"");
				appendStrW(data,Int2Wstr(outputWidth,tmp,10));
				appendStrW(data,L"\" biheight=\"");
				appendStrW(data,Int2Wstr(outputHeight,tmp,10));
				appendStrW(data,L"\"/>");
				appendStrW(data,L"</videoinfoheader>");
				appendStrW(data,L"</wmmediatype>");
				appendStrW(data,L"</streamconfig>");
				appendStrW(data,L"</profile>");
			}
	
			_RPT1(_CRT_WARN,"data:  %S",data);
	
			if (FAILED(hr = pWMProfileManager->LoadProfileByData(data,&pProfile))) mexErrMsgTxt("WMProfileManager can not load prx profile.  If you aren't loading a custom prx file, this most likely means that your video size isn't a multiple of 2.  Make sure both the width and height are even.");
			else _RPT0(_CRT_WARN,"LoadProfileByData succeeded!!!");
			
			delete [] data;
	
			pWMProfileManager->CreateEmptyProfile(WMT_VER_9_0,&pProfile2);
	
			DWORD count = 0;
			pProfile->GetStreamCount(&count);
			int streamNr = 1;
		
			for (int i=0; i<count; i++)
			{
				IWMStreamConfig* pConfig;
				pProfile->GetStream(i,&pConfig);
				
				GUID streamtype;
				pConfig->GetStreamType(&streamtype);
	
				if (streamtype == WMMEDIATYPE_Video)
				{
					for (int j=0; j<nrVideo; j++)
					{
						pConfig->SetStreamNumber(streamNr++);
						pProfile2->AddStream(pConfig);
					}
				} else if (streamtype == WMMEDIATYPE_Audio) {
					for (int j=0; j<nrAudio; j++)
					{
						pConfig->SetStreamNumber(streamNr++);
						pProfile2->AddStream(pConfig);
					}
				}
				pConfig->Release();
			}
	
			if (FAILED(hr = pConfigAsfWriter->ConfigureFilterUsingProfile(pProfile2))) mexErrMsgTxt("ConfigAsfWriter can not use prx profile.");
			else _RPT0(_CRT_WARN,"ConfigureFilterUsingProfile succeeded!!!");
			
			pConfigAsfWriter->Release();
			pProfile->Release();
			pProfile2->Release();
			pWMProfileManager->Release();
		}

		IFileSinkFilter *pSink= NULL;
		pWriter->QueryInterface(IID_IFileSinkFilter, (void**)&pSink);
		pSink->SetFileName(uniFilename, NULL);
		pSink->Release();
	
		// connect mmfilt to Writer/Mux+Writer
		IPin* pin = NULL;
		IEnumPins* pinList;
	
		//get the output pins so that we can render them...
		if (SUCCEEDED(mmfilt->EnumPins(&pinList)))
		{
			pinList->Reset();
			while (pinList->Next(1, &pin, NULL) == S_OK)
			{
				_RPT0(_CRT_WARN,"trying to insert video encoder\n");
				if (nrVideo==0 || FAILED(insertEncoder(pGraphBuilder, VideoEncoder, pMux, pin)))
				{
					_RPT0(_CRT_WARN,"trying to insert audio encoder\n");
					if (nrAudio==0 || FAILED(hr = insertEncoder(pGraphBuilder, AudioEncoder, pMux, pin)))
					{
						mexErrMsgTxt("Can't connect output pin as either Audio or Video");
					}
					_RPT2(_CRT_WARN,"hr %d %d\n",hr,FAILED(hr));
				}
				//pin->Release(); this is done by insertEncoder
			}
			pinList->Release();
		}

		IMediaFilter* pMediaFilter;
		pGraphBuilder->QueryInterface(IID_IMediaFilter, (void **)&pMediaFilter);
		//turn off the clock, so that it will run as fast as possible
		pMediaFilter->SetSyncSource(NULL);
		pMediaFilter->Release();
	}

	if (buildGraph)
	{
		pGraphBuilder->QueryInterface(IID_IMediaEvent, (void **)&pMediaEventEx);
		pGraphBuilder->QueryInterface(IID_IMediaControl, (void **)&pMediaControl);

		_RPT0(_CRT_WARN,"Trying to Run\n");
		if (FAILED(hr = pMediaControl->Run())) mexErrMsgTxt(message(hr));
		_RPT0(_CRT_WARN,"Running...\n");
	} else {
		//wake up the sleeping sources...
		for (int i=0;i<mmfilt->audioSources.size();i++) SetEvent(mmfilt->audioSources.at(i)->sleepEvent);
		for (int i=0;i<mmfilt->videoSources.size();i++) SetEvent(mmfilt->videoSources.at(i)->sleepEvent);
	}

	Sleep(100);

	long evCode = 0;
	int notDone = 1;
	while (pMediaEventEx->WaitForCompletion(1000, &evCode) == E_ABORT && notDone)
	{
		notDone = 0;
		for (int i=0;i<mmfilt->audioSources.size();i++) notDone |= mmfilt->audioSources.at(i)->done==0;
		for (int i=0;i<mmfilt->videoSources.size();i++) notDone |= mmfilt->videoSources.at(i)->done==0;
		_RPT1(_CRT_WARN,"notDone %d\n",notDone);
	}

	if (buildGraph)
	{
		if (pAviMux) pAviMux->Release();
		pWriter->Release();
	}
	if (stopGraph)
	{
		pMediaEventEx->Release();

		// make sure everything has really stopped before returning.
		if (FAILED(hr = pMediaControl->Stop())) mexErrMsgTxt(message(hr));

		pMediaControl->Release();

		cleanUp(&pGraphBuilder);
	}
	
//	delete mmfilt; // this is automatically deleted???
}

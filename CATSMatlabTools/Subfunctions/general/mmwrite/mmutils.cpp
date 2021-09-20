TCHAR str[200];

TCHAR* message(HRESULT hr)
{
	if (hr == S_OK)
	{
		return "";
	} else {
		if (AMGetErrorText(hr,str,200) != 0) return str;
		return "Unknown error";
	}
}

vector<void*> newlist;

void* newadd(void* item)
{
	newlist.add(item);
	return item;
}

void cleanUp(IGraphBuilder** pGraphBuilder)
{
	IMediaControl* pMediaControl;

	vector<IBaseFilter*> filts;
	IEnumFilters* filterList;
	IBaseFilter* filt;

	_RPT0(_CRT_WARN,"Releasing... \n");
	if (SUCCEEDED((*pGraphBuilder)->EnumFilters(&filterList)))
	{
		filterList->Reset();
		while (filterList->Next(1, &filt, NULL) == S_OK) filts.add(filt);
		filterList->Release();
	}
	
	for (int i=0;i<filts.size();i++)
	{
#ifdef _DEBUG
		FILTER_INFO info;
		filts.at(i)->QueryFilterInfo(&info);
		char str[100];
		WideCharToMultiByte( CP_ACP, 0, info.achName, -1, str, 100, NULL, NULL );
		_RPT1(_CRT_WARN,"Releasing: %s\n",str);
#endif
		(*pGraphBuilder)->RemoveFilter(filts.at(i));
		filts.at(i)->Release();
	}
	(*pGraphBuilder)->Release();
	(*pGraphBuilder) = NULL;

	for (int i=0;i<newlist.size();i++) if (newlist.at(i)) delete [] newlist.at(i);
	newlist.clear();

	_CrtDumpMemoryLeaks();
}

void MyFreeMediaType(AM_MEDIA_TYPE& mt)
{
	if (mt.cbFormat != 0)
	{
		CoTaskMemFree((PVOID)mt.pbFormat);
		mt.cbFormat = 0;
		mt.pbFormat = NULL;
	}
	if (mt.pUnk != NULL)
	{
		// Unecessary because pUnk should not be used, but safest.
		mt.pUnk->Release();
		mt.pUnk = NULL;
	}
}

PIN_INFO getPinInfo(IPin* pin)
{
	PIN_INFO	info = {0};

	if (pin)
	{
		if (SUCCEEDED(pin->QueryPinInfo(&info)))
		{
			info.pFilter->Release();
		}
	}

	return info;
}

IPin* getInputPin(IBaseFilter* filt)
{
	IPin* pin = NULL;
	IEnumPins* pinList;

	if (!filt) return NULL;

	//get the input
	if (SUCCEEDED(filt->EnumPins(&pinList)))
	{
		pinList->Reset();
		while (pinList->Next(1, &pin, NULL) == S_OK && getPinInfo(pin).dir != PINDIR_INPUT) pin->Release();
		pinList->Release();

		if (getPinInfo(pin).dir != PINDIR_INPUT) return NULL;
	}

	return pin;
}

IPin* getOutputPin(IBaseFilter* filt)
{
	IPin* pin = NULL;
	IEnumPins* pinList;

	if (!filt) return NULL;

	//get the output
	if (SUCCEEDED(filt->EnumPins(&pinList)))
	{
		pinList->Reset();
		while (pinList->Next(1, &pin, NULL) == S_OK && getPinInfo(pin).dir != PINDIR_OUTPUT) pin->Release();
		pinList->Release();

		if (getPinInfo(pin).dir != PINDIR_OUTPUT) return NULL;
	}

	return pin;
}

bool isRenderer(IBaseFilter* filt)
{
	if (!filt) return false;

	IEnumPins*	pinList;
	int nrOutput = 0;
	int nrInput = 0;
	IPin*		pin = NULL;

	if (FAILED(filt->EnumPins(&pinList))) return false;
	pinList->Reset();
	while (pinList->Next(1, &pin, NULL) == S_OK)
	{
		if (getPinInfo(pin).dir == PINDIR_OUTPUT) nrOutput++;
		else nrInput++;
		pin->Release();
	}
	pinList->Release();

	#ifdef _DEBUG
		FILTER_INFO info;
		filt->QueryFilterInfo(&info);
		char str[100];
		WideCharToMultiByte( CP_ACP, 0, info.achName, -1, str, 100, NULL, NULL );
		_RPT0(_CRT_WARN,str);
		_RPT2(_CRT_WARN," %d %d\n", nrOutput, nrInput);
	#endif

	return nrOutput == 0 && nrInput == 1;  // the only filters that have no outputs are renderers
}


// impement our own versions of the functions which require the MSVCR dlls...
LPWSTR Int2Wstr(int i, LPWSTR wstr, size_t len)
{
	TCHAR tmp[32];
	for (int j=0;j<len;j++)
	{
		tmp[j] = (i % 10) + L'0';
		i = i / 10;
	}
	int c = 0;
	for (int j=len-1;j>=0;j--) if (c > 0 || tmp[j] > L'0') wstr[c++] = tmp[j];
	// null termination
	if (c<len) wstr[c] = 0;
	else wstr[len-1] = 0;
	
	return wstr;
}


void appendStrW(LPWSTR dest, LPCWSTR src)
{
	int len = 0,i=0;
	while (dest[len++] != 0);
	len--;
	while (src[i++] != 0) dest[len++] = src[i-1];
	dest[len] = 0; //null
}
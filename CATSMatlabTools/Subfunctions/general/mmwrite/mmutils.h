// make our own vector class because we are trying not to link to the MS dlls
template<class T> class vector
{
	public:
		vector()
		{
			datavec = NULL;
			datavecSize = 0;
			nr = 0;
			datavecSize = resize(128);
		}

		~vector()
		{
			if (datavec) free(datavec);
		}

		T at(unsigned int i)
		{
			if (i >= nr) return NULL;
			return datavec[i];
		}

		int size()
		{
			return nr;
		}

		void assign(unsigned int i, T data)
		{
			if (i >= nr) return;
			datavec[i] = data;
		}

		void add(T data)
		{
			nr++;
			if (nr > datavecSize)
			{
				datavecSize = resize(datavecSize*2);
				if (nr > datavecSize)
				{
					//we've failed, but can't return an error here...
					nr--;
					return;
				}
			}
			datavec[nr-1] = data;
		}

		void clear()
		{
			nr = 0;
			datavecSize = resize(128);
		}

	private:
		T* datavec;
		unsigned int nr;
		unsigned int datavecSize;

		int resize(int newsize)
		{
			void* tmp = new T[newsize];
			if (tmp)
			{
				memcpy(tmp,datavec,nr*sizeof(T));
				datavec = (T*)tmp;
				return newsize;
			}
			return datavecSize;
		}
};

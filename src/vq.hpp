#pragma once

namespace vq
{
	// Calculate distance between two nodes
	template <class T>
	double distance(T* a, T* b, int aNodeSize)
	{
		double d = 0;
		for (int i = 0; i < aNodeSize; i++)
			d += (static_cast<double>(a[i]) - static_cast<double>(b[i])) * (static_cast<double>(a[i]) - static_cast<double>(b[i]));
		return sqrt(d);
	}

	// Select initial bins (center nodes)
	template <class T>
	void select(T* aNode, int aNodeSize, int aNodeCount, T* aBin, int bins)
	{
		return;
	}

	// Average bin to get to a new center
	template <class T>
	void average(T* aNode, int aNodeSize, int* index, int indices, T* out)
	{
		if (!indices)
		{
			for (int i = 0; i < aNodeSize; i++)
				out[i] = 0;
			return;
		}

		for (int i = 0; i < aNodeSize; i++)
		{
			long long sum = 0;
			for (int j = 0; j < indices; j++)
			{
				sum += static_cast<long long>(aNode[index[j] * aNodeSize + i]);
			}
			sum /= indices;
			out[i] = (T)(sum);
		}
	}

	// Nudge node to split a bin
	template <class T>
	void nudge(T* aNode, int aNodeSize)
	{
		int i = rand() % aNodeSize;
		if (aNode[i] >= 1)
		{
			if (rand() & 1)
				aNode[i]--;
			else
				aNode[i]++;
		}
		else
	 	    aNode[i]++;
	}

	template <class T, int bins, int errorTreshold = 0, int maxiter = 1000000>
	void reduce(T* aNode, int aNodeSize, int aNodeCount, T* aBin)
	{
		int* bin;
		int binc[bins];
		for (auto& x : binc) x = 0;
		bin = new int[bins * aNodeCount];

		// 1. initial selection

		select(aNode, aNodeSize, aNodeCount, aBin, bins);

		double lastError = 1e100;
		double error = 1e100;
		int iter = 0;
		int force = 0;

		do 
		{
			force = 0;
			iter++;

			// 2. move to closest group

			for (auto& x : binc) x = 0;

			for (int ni = 0; ni < aNodeCount; ni++)
			{
				int closest = 0;
				double mindist = distance(aNode + ni * aNodeSize, aBin, aNodeSize);
				for (int bi = 0; bi < bins; bi++)
				{
					double d = distance(
						aNode + ni * aNodeSize, 
						aBin + bi * aNodeSize, 
						aNodeSize);
					if (d < mindist)
					{
						mindist = d;
						closest = bi;
					}
				}
				bin[closest * aNodeCount + binc[closest]] = ni;
				binc[closest]++;
			}
						
			// 3. calculate new centers

			for (int i = 0; i < bins; i++)
				if (binc[i])
					average(
						aNode, 
						aNodeSize, 
						bin + i * aNodeCount, 
						binc[i], 
						aBin + i * aNodeSize);

			// 4. measure error

			lastError = error;
			error = 0;
			for (int i = 0; i < bins; i++)
			{
				for (int j = 0; j < binc[i]; j++)
					error += distance(
						aBin + i * aNodeSize, 
						aNode + bin[i * aNodeCount + j] * aNodeSize, 
						aNodeSize);
			}

			printf("%c\r", "\\-/|"[iter % 4]);
			
			int emptybin = -1;
			for (int i = 0; i < bins; i++)
			{
				if (binc[i] == 0)
					emptybin = i;
			}

			if (emptybin != -1)
			{
				force = aNodeCount > bins;
				// there's an empty bin (or more, but let's fill one at a time)
				// - find biggest (not the most populated) group
				int biggest = 0;
				double bigsize = 0;
				for (int i = 0; i < bins; i++)
					for (int j = 0; j < binc[i]; j++)
					{
						double l = distance(aBin + i * aNodeSize, aNode + bin[i * aNodeCount + j] * aNodeSize, aNodeSize);
						if (l > bigsize)
						{
							biggest = i;
							bigsize = l;
						}
					}
				// - compete with its center
				for (int i = 0; i < aNodeSize; i++)
					aBin[emptybin * aNodeSize + i] = aBin[biggest * aNodeSize + i];
				nudge(aBin + emptybin * aNodeSize, aNodeSize);
			}
			
			// 5. repeat

		} 
		while (iter < maxiter && ((force && error) || (error - lastError) < errorTreshold));
		delete[] bin;
	}
}

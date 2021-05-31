#ifndef SIGNALSMITH_DSP_COMMON_H
#define SIGNALSMITH_DSP_COMMON_H

#include <cmath>

namespace signalsmith {
	namespace windows {
		/// Used for the Kaiser window
		inline static double bessel0(double x) {
			constexpr double significanceLimit = 1e-6;
			double result = 0;
			double term = 1;
			double m = 0;
			while (term  > result*significanceLimit) {
				result += term;
				++m;
				term *= (x*x) / (4*m*m);
			}

			return result;
		}
		/// Fills a container with a Kaiser window, with shape parameter `beta = alpha*pi`
		template<typename Data>
		void fillKaiser(Data &data, int size, double beta) {
			double scaling = 1/bessel0(beta);
			double invSize = 1.0/size;
			for (int i = 0; i < size; ++i) {
				double r = (2*i + 1)*invSize - 1;
				double arg = std::sqrt(1 - r*r);
				data[i] = bessel0(beta*arg)*scaling;
			}
		}
		/** Fills a container with a Kaiser window with a specified bandwidth (as a factor of 1/window-length)
		
		Bandwidth is measured as the the main lobe width (between nulls), unless `approximateOptimal` is used, in which case the boundary is placed _slightly_ inside the main-lobe, following an approximate compromise between total/peak energy ratios on either side of the boundary.*/
		template<typename Data>
		void fillKaiserBandwidth(Data &data, int size, double bandwidth, bool approximateOptimal) {
			if (approximateOptimal) {
				// Heuristic based on numerical search
				bandwidth = bandwidth + 2/(bandwidth*bandwidth);
			}
			double alpha = std::sqrt(bandwidth*bandwidth*0.25 - 1);
			fillKaiser(data, size, alpha*M_PI);
		}
		/// Roughly optimal Kaiser for STFT analysis (forced to perfect reconstruction)
		template<typename Data>
		void fillKaiserStft(Data &data, int size, int stride) {
			double overlap = size/(double)stride;
			fillKaiserBandwidth(data, size, overlap, true);
			
			// Force perfect reconstruction
			for (int i = 0; i < stride; ++i) {
				double sum2 = 0;
				for (int index = i; index < size; index += stride) {
					sum2 += data[index]*data[index];
				}
				double factor = std::sqrt(1/sum2);
				for (int index = i; index < size; index += stride) {
					data[index] *= factor;
				}
			}
		}
	}

} // signalsmith::
#endif // include guard

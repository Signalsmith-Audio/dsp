#ifndef SIGNALSMITH_DSP_WINDOWS_H
#define SIGNALSMITH_DSP_WINDOWS_H

#include "./common.h"

#include <cmath>

namespace signalsmith {
namespace windows {
	/**	@defgroup Windows Window functions
		@brief Windows for spectral analysis
		
		These are generally double-precision, because they are mostly calculated during setup/reconfiguring, not real-time code.
		
		@{
		@file
	*/
	
	/** @brief The Kaiser window (almost) maximises the energy in the main-lobe compared to the side-lobes.
		These can be specified by the shape-parameter (beta) or by the main lobe's bandwidth (with an optional heuristic: see: `.withBandwidth()`):
		\diagram{kaiser-windows.svg,You can see that the main lobe matches the specified bandwidth.}
	*/
	class Kaiser {
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
		
		double beta;
		double invB0;
	public:
		/// Set up a Kaiser window with a given shape.  `beta` is `pi*alpha` (since there is ambiguity about shape parameters)
		Kaiser(double beta) : beta(beta), invB0(1/bessel0(beta)) {}
		
		/** @brief Returns a Kaiser window where the main lobe has the specified bandwidth (as a factor of 1/window-length).
		If `approximateOptimal` is enabled, the main lobe width is _slightly_ wider, following an approximate compromise between total/peak energy ratios on either side of the boundary.
		\diagram{kaiser-windows-heuristic.svg,Compare this to the diagram at the top\, where the main-lobe ends exactly at the boundary.}
		*/
		static Kaiser withBandwidth(double bandwidth, bool approximateOptimal=false) {
			if (approximateOptimal) {
				// Heuristic based on numerical search
				bandwidth = bandwidth + 2/(bandwidth*bandwidth); // mostly improves peak/avg, but gives slightly higher peaks for small bandwidth (in exchange for better average case)
			}
			double alpha = std::sqrt(bandwidth*bandwidth*0.25 - 1);
			return Kaiser(alpha*M_PI);
		}

		/// Fills an arbitrary container with a Kaiser window
		template<typename Data>
		void fill(Data &data, int size) const {
			double invSize = 1.0/size;
			for (int i = 0; i < size; ++i) {
				double r = (2*i + 1)*invSize - 1;
				double arg = std::sqrt(1 - r*r);
				data[i] = bessel0(beta*arg)*invB0;
			}
		}
		
	};

	/**@brief Forces STFT perfect-reconstruction (WOLA) on an existing window by direct calculation.
	For example, here are perfect-reconstruction versions of the approximately-optimal @ref Kaiser windows:
	\diagram{kaiser-windows-heuristic-pr.svg,Note the lower overall energy\, and the pointy top for 2x bandwidth. Spectral performance is about the same, though.}
	*/
	template<typename Data>
	void forcePerfectReconstruction(Data &data, int size, int stride) {
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

/** @} */
}} // signalsmith::windows
#endif // include guard

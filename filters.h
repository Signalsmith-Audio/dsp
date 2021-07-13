#ifndef SIGNALSMITH_DSP_SPECTRAL_H
#define SIGNALSMITH_DSP_SPECTRAL_H

#include "./common.h"

#include "./fft.h"
#include "./delay.h"

#include <cmath>

namespace signalsmith {
namespace filters {
	/**	@defgroup Filters Basic filters
		@brief Classes for some common filter types
		
		@{
		@file
	*/
	
	/** A standard biquad - not guaranteed to be stable if modulated at audio rate */
	template<typename Sample, bool correctBandwidth=false>
	class BiquadStatic {
		Sample a1 = 0, a2 = 0, b0 = 1, b1 = 0, b2 = 0;
		Sample x1 = 0, x2 = 0, y1 = 0, y2 = 0;
		
		// Straight from the cookbook: https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
		enum class Type {highpass, lowpass, highShelf, lowShelf, bandpass, bandStop};
		SIGNALSMITH_INLINE void configure(Type type, double scaledFreq, double octaves, double db=-6) {
			if (scaledFreq >= 0.499) scaledFreq = 0.499;
			double w0 = 2*M_PI*scaledFreq;
			double cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
			double alpha = sin_w0*std::sinh(std::log(2)*0.5*octaves*(correctBandwidth ? w0/sin_w0 : 1));
			double A = std::pow(10, db/40), sqrtA2alpha = std::sqrt(A)*alpha;

			double a0;
			if (type == Type::highpass) {
				b1 = -1 - cos_w0;
				b0 = b2 = (1 + cos_w0)*0.5;
				a0 = 1 + alpha;
				a1 = -2*cos_w0;
				a2 = 1 - alpha;
			} else if (type == Type::lowpass) {
				b1 = 1 - cos_w0;
				b0 = b2 = b1*0.5;
				a0 = 1 + alpha;
				a1 = -2*cos_w0;
				a2 = 1 - alpha;
			} else if (type == Type::highShelf) {
				b0 = A*((A+1)+(A-1)*cos_w0+sqrtA2alpha);
				b2 = A*((A+1)+(A-1)*cos_w0-sqrtA2alpha);
				b1 = -2*A*((A-1)+(A+1)*cos_w0);
				a0 = (A+1)-(A-1)*cos_w0+sqrtA2alpha;
				a2 = (A+1)-(A-1)*cos_w0-sqrtA2alpha;
				a1 = 2*((A-1)-(A+1)*cos_w0);
			} else if (type == Type::lowShelf) {
				b0 = A*((A+1)-(A-1)*cos_w0+sqrtA2alpha);
				b2 = A*((A+1)-(A-1)*cos_w0-sqrtA2alpha);
				b1 = 2*A*((A-1)-(A+1)*cos_w0);
				a0 = (A+1)+(A-1)*cos_w0+sqrtA2alpha;
				a2 = (A+1)+(A-1)*cos_w0-sqrtA2alpha;
				a1 = -2*((A-1)+(A+1)*cos_w0);
			} else if (type == Type::bandpass) {
				b0 = alpha;
				b1 = 0;
				b2 = -alpha;
				a0 = 1 + alpha;
				a1 = -2*cos_w0;
				a2 = 1 - alpha;
			} else if (type == Type::bandStop) {
				b0 = 1;
				b1 = -2*cos_w0;
				b2 = 1;
				a0 = 1 + alpha;
				a1 = b1;
				a2 = 1 - alpha;
			} else {
				// reset to neutral
				a1 = a2 = b1 = b2 = 0;
				a0 = b0 = 1;
			}
			double invA0 = 1/a0;
			b0 *= invA0;
			b1 *= invA0;
			b2 *= invA0;
			a1 *= invA0;
			a2 *= invA0;
		}
	public:

		Sample operator ()(Sample x0) {
			Sample y0 = x0*b0 + x1*b1 + x2*b2 - y1*a1 - y2*a2;
			y2 = y1;
			y1 = y0;
			x2 = x1;
			x1 = x0;
			return y0;
		}
		
		void reset() {
			x1 = x2 = y1 = y2 = 0;
		}

		void highpass(double scaledFreq, double octaves=2) {
			configure(Type::highpass, scaledFreq, octaves);
		}
		void lowpass(double scaledFreq, double octaves=2) {
			configure(Type::lowpass, scaledFreq, octaves);
		}
		void highShelf(double scaledFreq, double octaves=2, double db=-6) {
			configure(Type::highShelf, scaledFreq, octaves, db);
		}
		void lowShelf(double scaledFreq, double octaves=2, double db=-6) {
			configure(Type::lowShelf, scaledFreq, octaves, db);
		}
		void bandpass(double scaledFreq, double octaves=1) {
			configure(Type::bandpass, scaledFreq, octaves);
		}
		void bandStop(double scaledFreq, double octaves=1) {
			configure(Type::bandStop, scaledFreq, octaves);
		}
	};

/** @} */
}} // signalsmith::filters::
#endif // include guard

#ifndef SIGNALSMITH_DSP_FILTERS_H
#define SIGNALSMITH_DSP_FILTERS_H

#include "./common.h"
#include "./perf.h"

#include <cmath>

namespace signalsmith {
namespace filters {
	/**	@defgroup Filters Basic filters
		@brief Classes for some common filter types
		
		@{
		@file
	*/
	
	enum class BiquadDesign {
		bilinear, /// Bilinear transform, adjusting for centre frequency but not bandwidth
 		cookbook, /// RBJ's "Audio EQ Cookbook".  Based on `bilinear`, adjusts bandwidth (for peak/notch/bandpass) by preserving the ratio between upper/lower boundaries.  This performs oddly near Nyquist.
		oneSided, /// Based on `bilinear`, but adjust bandwidth by preserving the lower boundary (and leaves the upper one loose).
		vicanek /// From Martin Vicanek's "Matched Second Order Digital Filters".  This takes the poles from the impulse-invariant approach, and then picks the zeros to create a better match.  This means that Nyquist is not 0dB for peak/notch (or -Inf for lowpass), but it is a decent match to the analogue prototype. https://vicanek.de/articles/BiquadFits.pdf
	};
	
	/** A standard biquad.

		This is not guaranteed to be stable if modulated at audio rate.
		
		The default highpass/lowpass bandwidth (1.9) produces a Butterworth filter when bandwidth-compensation is disabled.
		
		Bandwidth compensation defaults to `true` for all filter types aside from highpass/lowpass.  When in "cookbook mode", it roughly preserves the ratio between upper/lower boundaries (-3dB point or half-gain).  Otherwise, it exactly preserves the lower boundary.*/
	template<typename Sample, bool cookbookBandwidth=false>
	class BiquadStatic {
		static constexpr BiquadDesign bwDesign = cookbookBandwidth ? BiquadDesign::cookbook : BiquadDesign::oneSided;
		Sample a1 = 0, a2 = 0, b0 = 1, b1 = 0, b2 = 0;
		Sample x1 = 0, x2 = 0, y1 = 0, y2 = 0;
		
		// Straight from the cookbook: https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
		enum class Type {highpass, lowpass, highShelf, lowShelf, bandpass, bandStop};

		SIGNALSMITH_INLINE void configure(Type type, double scaledFreq, double octaves, double sqrtGain, BiquadDesign design) {
			scaledFreq = std::max(0.0001, std::min(0.4999, scaledFreq));
			double w0 = 2*M_PI*scaledFreq;
			if (design == BiquadDesign::vicanek) {
				double Q = 0.5/std::sinh(std::log(2)*0.5*octaves);
				double q = 1/(2*Q);
				double expmqw = std::exp(-q*w0);
				if (q <= 1) {
					a1 = -2*expmqw*std::cos(std::sqrt(1 - q*q)*w0);
				} else {
					a1 = -2*expmqw*std::cosh(std::sqrt(q*q - 1)*w0);
				}
				a2 = expmqw*expmqw;
				double sinpd2 = std::sin(w0/2);
				double p0 = 1 - sinpd2*sinpd2, p1 = sinpd2*sinpd2, p2 = 4*p0*p1;
				double A0 = 1 + a1 + a2, A1 = 1 - a1 + a2, A2 = -4*a2;
				A0 *= A0;
				A1 *= A1;
				if (type == Type::lowpass) {
					double R1 = (A0*p0 + A1*p1 + A2*p2)*Q*Q;
					double B0 = A0, B1 = (R1 - B0*p0)/p1;
					b0 = 0.5*(std::sqrt(B0) + std::sqrt(B1));
					b1 = std::sqrt(B0) - b0;
					b2 = 0;
					return;
				} else if (type == Type::highpass) {
					b2 = b0 = std::sqrt(A0*p0 + A1*p1 + A2*p2)*Q/(4*p1);
					b1 = -2*b0;
					return;
				}
			}
			double cos_w0 = std::cos(w0), sin_w0 = std::sin(w0);
			if (design != BiquadDesign::bilinear) {
				if (design == BiquadDesign::cookbook) {
					// Approximately preserves bandwidth between halfway points
					octaves *= w0/sin_w0;
				} else {
					// Bilinear warps frequencies like tan(pi*x)/pi, so this places the lower boundary in the correct place
					double lowerRatio = std::pow(2, octaves*-0.5);
					octaves = 2*std::log2(std::tan(M_PI*scaledFreq)/std::tan(M_PI*scaledFreq*lowerRatio));
				}
			}
			double alpha = sin_w0*std::sinh(std::log(2)*0.5*octaves);
			double A = sqrtGain, sqrtA2alpha = std::sqrt(A)*alpha;

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

		void lowpass(double scaledFreq, BiquadDesign design=BiquadDesign::bilinear) {
			return lowpass(scaledFreq, 1.9, design);
		}
		void lowpass(double scaledFreq, double octaves, BiquadDesign design=BiquadDesign::bilinear) {
			configure(Type::lowpass, scaledFreq, octaves, 0, design);
		}
		void highpass(double scaledFreq, BiquadDesign design=BiquadDesign::bilinear) {
			return highpass(scaledFreq, 1.9, design);
		}
		void highpass(double scaledFreq, double octaves, BiquadDesign design=BiquadDesign::bilinear) {
			configure(Type::highpass, scaledFreq, octaves, 0, design);
		}

		// Old API
		void highpass(double scaledFreq, double octaves, bool correctBandwidth) {
			configure(Type::highpass, scaledFreq, octaves, 0, correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void lowpass(double scaledFreq, double octaves, bool correctBandwidth) {
			configure(Type::lowpass, scaledFreq, octaves, 0, correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void bandpass(double scaledFreq, double octaves=1, bool correctBandwidth=true) {
			configure(Type::bandpass, scaledFreq, octaves, 0, correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void highShelf(double scaledFreq, double gain, double octaves=2, bool correctBandwidth=true) {
			configure(Type::highShelf, scaledFreq, octaves, std::sqrt(gain), correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void lowShelf(double scaledFreq, double gain, double octaves=2, bool correctBandwidth=true) {
			configure(Type::lowShelf, scaledFreq, octaves, std::sqrt(gain), correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void highShelfDb(double scaledFreq, double db, double octaves=2, bool correctBandwidth=true) {
			double sqrtGain = std::pow(10, db*0.025);
			configure(Type::highShelf, scaledFreq, octaves, sqrtGain, correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void lowShelfDb(double scaledFreq, double db, double octaves=2, bool correctBandwidth=true) {
			double sqrtGain = std::pow(10, db*0.025);
			configure(Type::lowShelf, scaledFreq, octaves, sqrtGain, correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
		void bandStop(double scaledFreq, double octaves=1, bool correctBandwidth=true) {
			configure(Type::bandStop, scaledFreq, octaves, 0, correctBandwidth ? bwDesign : BiquadDesign::bilinear);
		}
	};

/** @} */
}} // signalsmith::filters::
#endif // include guard

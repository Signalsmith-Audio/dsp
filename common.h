#ifndef SIGNALSMITH_DSP_COMMON_H
#define SIGNALSMITH_DSP_COMMON_H

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

namespace signalsmith {
	/**	@defgroup Common Common
		@brief helper classes used by the rest of the library
		
		@{
		@file
	*/
//	
//	/// Helpful base-class for anything which could be either fixed-size or variable size (`-1`)
//	template<int size=-1>
//	class Size {
//	public:
//		Size() {}
//
//		static constexpr int size() {
//			return size;
//		}
//	}
//
//	/// Variable-size specialisation
//	template<>
//	class Size<-1> {
//		int _size;
//	public:
//		Size(int size=0) : _size(size) {}
//
//		int size() const {
//			return _size;
//		}
//	}

/** @} */
} // signalsmith::
#endif // include guard

# Signalsmith Audio's DSP Library

A C++11 header-only library, providing classes/templates for (mostly audio) signal-processing tasks.

```
git clone https://signalsmith-audio.co.uk/code/dsp.git
```

Just include the header file(s) you need, and an optional version-check:

```
#include "dsp/delay.h"
SIGNALSMITH_DSP_VERSION_CHECK(1, 0, 0)
```

### Documentation

Docs are on the web: https://signalsmith-audio.co.uk/code/dsp/

### Development / contributing

Tests (and source-scripts for the above docs) are available in a separate repo:

```
git clone https://signalsmith-audio.co.uk/code/dsp-doc.git
```

The goal (where possible) is to measure/test the actual audio characteristics of the tools (e.g. frequency responses and aliasing levels).

### License

This code is [MIT licensed](LICENSE.txt).  If you'd prefer something else, get in touch.

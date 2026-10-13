PAX Scintillator Processing
---

In the PAX experiment we monitored beam intensity with a pair of square scintillators (30 cm x 30 cm) whose signal was digitized by a LeCroy Waverunner 8104.

I have here a quick pile of code I used to process the scintillator signals.
Final analysis needed some additional data sorting and matchting; there was a second roll-over problem and the clocks between the TES DAQ were not shared (TES system uses a 10 MHz Rb standard) until part way through the run.

`processs_scint.py`

This has the functions, including plotting tools, I used.
It requires the proprietary LecroyParser library. 
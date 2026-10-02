We needed an applet to watch the temperature of a cryo system.
There already was a tool to dump the temp value, time and the % of max heater output current to a csv live.

The bones of the animation were put together by Goncalo Baptista but I stole it and reworked some of the vis so now it's mine. 


plot_live.py
--

This plots (as written) the last 10 mins of temp evolution with a rolling time window updated every 2 seconds.
If I've done it right, this does not hold anything in the cache and has at least 1 night of over-night operations confirmed.


Example output: 
`python plot_live.py ADRLog_example.txt`

![Example](./Figure_1.png)


I'll note the computer this runs on has pandas 3.0.0 and my personal computer has 3.0.1. To get this to run with 3.0.1 you will probably need to add a delimiter flag to catch the `\t`s.  

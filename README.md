# POSSUM_channel_flagger
A channel-flagging script for the POSSUM pipeline.

To identify bad channels, two quantities are computed per channel:
* The most extreme value present (the maximum of the absolute value over all pixels in a channel).
* The channel noise, computed as the median absolute deviation from the median (MADFM), corrected to 
be equivalent to a Gaussian sigma.

Using these two quantities, two sets of tests are performed per-channel:
* excessively large values: a channel with an incredibly large extreme value or noise, in absolute
 terms, is considered bad. The user can supply the thresholds for what is considered incredibly
large (with separate thesholds for extreme value and noise).

* values are much larger or smaller than median of neighbouring channels.

To test against neighbouring channels, the local medians of extreme value and noise are computed.
The user sets the number of channels, N, and the channels between -N/2 to +N/2 of the current channel
are used to compute a local median. At the edges of the array, this gets clipped, so the edge channels
are considered to have fewer neighbours.

A threshold multiple (separately for the extreme and noise) is provided by the user.
The exact test is whether the channel value (extreme or noise) is larger than the local median
times that multiple (e.g., if the noise multiple is 5 then the test is whether the channel's noise
is greater than 5x the local median noise); channels above the threshold are flagged.
The same test is also done in reverse: if the channel value is below the local value *divided* by
the multiple (e.g., below 1/5th the local median) it is also flagged. The net result is that good channels
must have a value between 1/M and M times the local median (i.e., equal upper and lower bounds in log-space).

All of these tests are run per-channel: large absolute value, and local value, for extrema and noise.
If one or more tests fail, the channel is flagged. These tests are run for all available Stokes parameters
(currently assuming I,Q,U are all present; V is optional). If a channel fails for one or more Stokes parameter,
it is flagged for all Stokes parameters.

The default behaviour of the script is to output a file with the flags (1/True = good channel, 
0/False = bad channel), along with columns containing the extemes and noise for all 4 Stokes.
To actually flag the cubes (set the bad channels to NaNs), the -f keyword must be used (or the relevant functions).




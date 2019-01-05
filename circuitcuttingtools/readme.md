# Introduction #

This is a set of command line driven (actually [twistcodewiki](https://bitbucket.org/goatchurch/twistcodewiki))
circuit cutting tools that depends on the [flatcam](http://flatcam.org/) code to parse Gerber and Excellon files, 
and which in turn depends on [shapely](http://toblerity.org/shapely/manual.html) for its 2D offsetting.

Since this is just 2D follow the dots CAM I don't need much UI (except for preview), but I do require some 
additional hackability on the way that the toolpaths come out.  If I do many boards with the same tooling I 
should be able to really automate it.

Features to add in when doing the next job:

* probing out the copper and then readproberelay.py to reset the z-heights to cut the groove on the warped plane [DONE]

* Automatic spot drilling the drillholes with the V-groove cutter while in the chuck [DONE]

* better choice of start positions on a contour so it's not right on a corner [partyly DONE]

* running along beyond the end-point to ensure isolation.  [DONE]
(I really need to run some cutting trials to prove this effect of the copper bridge when the start and end coincides -- 
and then file this as an issue to flatcam as a necessary feature.  blog any pictures)



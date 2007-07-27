There have been many scripts. All in repository history. Now, we keep here
the key ones:

tgd.sequential.R for the TGD in its sequential form, with original code
(will use gd1.R).

tgd.parallel.new.R. Run as 
nohup R --slave --no-save --no-restore --args 1 <tgd.parallel.new.R
>tgd.parallel.new.1.Rout

modifying the number for the number of slaves


and 

fcms.paral.R which is also run as above.


jclbasic.ho: jclbasic.h
jclbasic.ho: /home/jstolfi/include/bool.ho

jclparse.ho: jclparse.h
jclparse.ho: jclbasic.ho

jclsort.ho: jclsort.h

jclimage.ho: jclimage.h
jclimage.ho: jclbasic.ho

jcldist.ho: jcldist.h
jcldist.ho: jclbasic.ho

jclerror.ho: jclerror.h

jcltree.ho: jcltree.h
jcltree.ho: jclbasic.ho

jclgeo.ho: jclgeo.h

jcloptions.ho: jcloptions.h
jcloptions.ho: jclbasic.ho
jcloptions.ho: /home/jstolfi/include/bool.ho

jclbasic.o: jclbasic.c
jclbasic.o: jclbasic.ho
jclbasic.o: jclerror.ho

jclparse.o: jclparse.c
jclparse.o: jclbasic.ho
jclparse.o: jclerror.ho

jclsort.o: jclsort.c

geocluster.o: geocluster.c
geocluster.o: /home/jstolfi/include/bool.ho
geocluster.o: /home/jstolfi/include/argparser.ho
geocluster.o: jclbasic.ho
geocluster.o: jcldist.ho
geocluster.o: jclerror.ho
geocluster.o: jclgeo.ho
geocluster.o: jclimage.ho
geocluster.o: jcloptions.ho
geocluster.o: jclparse.ho
geocluster.o: jclsort.ho
geocluster.o: jcltree.ho

jcldist.o: jcldist.c

jclerror.o: jclerror.c
jclerror.o: jclerror.ho
jclerror.o: jclbasic.ho

jclimage.o: jclimage.c
jclimage.o: jclimage.ho
jclimage.o: jclbasic.ho
jclimage.o: jclerror.ho

jcltree.o: jcltree.c
jcltree.o: jcltree.ho
jcltree.o: jclbasic.ho
jcltree.o: jclerror.ho

jclgeo.o: jclgeo.c
jclgeo.o: jclbasic.ho

jcloptions.o: jcloptions.c
jcloptions.o: jcloptions.ho
jcloptions.o: jclbasic.ho
jcloptions.o: jclerror.ho


<TeXmacs|1.0.7.19>

<style|generic>

<\body>
  <doc-data|<doc-title|DoGPack Documentation File>>

  <block|<tformat|<table|<row|<cell|<with|font-series|bold|Directory>:
  <with|font-family|tt|apps/2d/advection/smooth_example>>>|<row|<cell|<with|font-series|bold|Equation:>
  <math|q<rsub|,t>+u*q<rsub|,x>+v*q<rsub|,y>>=0>>>>>

  <section|Mathematics>

  Consider the space-time domain <math|<around*|[|0,T|]>\<times\>\<Omega\> >,
  where <math|\<Omega\>\<subseteq\>\<bbb-R\><rsup|2>\<nocomma\>>, and let
  <math|t\<in\><around*|[|0,T|]>\<nosymbol\>\<nosymbol\>> and
  <math|<around*|(|x,y|)>\<in\>\<Omega\>>. The dependent variable is the
  concentration of a passive tracer, denoted by
  <math|q:\<bbb-R\>\<times\>\<Omega\>\<rightarrow\>\<bbb-R\>>, which is
  advected by the velocity field <math|<around*|(|u,v|)>\<nocomma\>>, where
  <math|><math|u:\<Omega\>\<rightarrow\>\<bbb-R\>> is the <math|x>-component
  of the velocity field and <math|v:\<Omega\>\<rightarrow\>\<bbb-R\>> is the
  <math|y>-component of the velocity field. The concentration then satisfies
  the following time-dependent advection equation:

  <\eqnarray*>
    <tformat|<table|<row|<cell|q<rsub|,t>+u*q<rsub|,x>+v*q<rsub|,y>=0.>|<cell|>|<cell|>>>>
  </eqnarray*>

  From the method of characteristics, we find that
  <math|q<around*|(|t,x,y|)>> is constant along the characteristics that are
  defined by the velocity field:

  <\eqnarray*>
    <tformat|<table|<row|<cell|q<around*|(|t,x,y|)>=q<around*|(|0,\<xi\>,\<eta\>|)>\<nocomma\>\<comma\>
    <space|1em>where<space|1em><frac|d\<xi\>|d
    t>=u<around*|(|\<xi\>,\<eta\>|)>,<space|1em><frac|d\<eta\>|d
    t>=v<around*|(|\<xi\>,\<eta\>|)>>|<cell|>|<cell|with<space|1em>\<xi\><around*|(|0|)>=x,<space|1em>\<eta\><around*|(|0|)>=y.>>>>
  </eqnarray*>

  The affect of boundary conditions can also be readily incorporated into the
  method of characteristics.

  <section|Editing Basic Application Files>

  Below we describe all the basic files that can be edited to change basi

  <\enumerate-numeric>
    <item>The initial conditions can be changed in the following file:

    <\session|shell|default>
      <\input|Shell] >
        open $DOGPACK/apps/2d/advection/smooth_example/QinitFunc.cpp
      </input>

      <\input|Shell] >
        \;
      </input>
    </session>

    <item>The velocity field can be changed in the following file:

    <\session|shell|default>
      <\input|Shell] >
        open $DOGPACK/apps/2d/advection/smooth_example/AuxFunc.cpp
      </input>

      <\input|Shell] >
        \;
      </input>
    </session>

    <item>The boundary conditions can be changed in the following file:

    <\session|shell|default>
      <\input|Shell] >
        open $DOGPACK/apps/2d/advection/smooth_example/SetBndValues.cpp
      </input>

      <\input|Shell] >
        \;
      </input>
    </session>

    <item>The flux function can be changed in the following file:

    <\session|xypic|default>
      <\session|shell|default>
        <\input|Shell] >
          open $DOGPACK/apps/2d/advection/smooth_example/FluxFunc.cpp
        </input>

        <\input|Shell] >
          \;
        </input>
      </session>
    </session>

    <item>The wave speeds needed to determine the time-step can be changed in
    the following file:

    <\session|shell|default>
      <\input|Shell] >
        open $DOGPACK/apps/2d/advection/smooth_example/SetWaveSpd.cpp
      </input>

      <\input|Shell] >
        \;
      </input>
    </session>

    <item>The projection onto the left and right eigenvectors can be done in
    the following two files (these are only needed if moment-limiters are
    used):

    <\session|shell|default>
      <\input|Shell] >
        open $DOGPACK/apps/2d/advection/smooth_example/ProjectLeftEig.cpp
      </input>

      <\input|Shell] >
        open $DOGPACK/apps/2d/advection/smooth_example/ProjectRightEig.cpp
      </input>

      <\input|Shell] >
        \;
      </input>
    </session>

    \;
  </enumerate-numeric>

  <section|Compilation>

  The example in this directory can be compiled as follows:

  <\session|shell|default>
    <\output>
      Shell session inside TeXmacs pid = 5559
    </output>

    <\input|Shell] >
      cd $DOGPACK/apps/2d/advection/smooth_example/
    </input>

    <\input|Shell] >
      make \<gtr\>/dev/null 2\<gtr\>&1
    </input>

    <\unfolded-io|Shell] >
      ls dog.exe
    <|unfolded-io>
      dog.exe
    </unfolded-io>

    <\input|Shell] >
      \;
    </input>
  </session>

  <section|Editing parameter file>

  The example parameters can be changed in the parameter file:

  <\session|shell|default>
    <\input|Shell] >
      open $DOGPACK/apps/2d/advection/smooth_example/parameters.ini
    </input>

    <\input|Shell] >
      \;
    </input>
  </session>

  <section|Running the example>

  To run the example simply type:

  <\session|shell|default>
    <\input|Shell] >
      cd $DOGPACK/apps/2d/advection/smooth_example/
    </input>

    <\input|Shell] >
      dog.exe \<gtr\>/dev/null 2\<gtr\>&1
    </input>

    <\input|Shell] >
      \;
    </input>
  </session>

  <section|Plotting in Python>

  To plot the solution results in Python type the following:

  <\session|python|default>
    <\input|Python] >
      import os
    </input>

    <\input|Python] >
      DOG = os.environ["DOGPACK"];
    </input>

    <\input|Python] >
      wdir = "".join((DOG,"/apps/2d/advection/smooth_example"));
    </input>

    <\input|Python] >
      os.chdir(wdir);
    </input>

    <\input|Python] >
      import plotdog2np as pd
    </input>

    <\unfolded-io|Python] >
      pd.plotdog2np(4,'output',1,1);
    <|unfolded-io>
      GridType = \ Cartesian

      \ \ \ \ \ \ \ \ points_per_dir = \ 4

      \ \ \ \ \ \ \ \ \ \ \ \ point_type = \ 1

      \ \ \ \ \ \ \ \ \ \ \ \ \ outputdir = \ output

      \ component_of_solution = \ 1

      \;

      \ Finished creating file for FRAME = \ 0 \ \ filename =
      \ figure0000.jpg

      \ Finished creating file for FRAME = \ 1 \ \ filename =
      \ figure0001.jpg

      \ Finished creating file for FRAME = \ 2 \ \ filename =
      \ figure0002.jpg

      \ Finished creating file for FRAME = \ 3 \ \ filename =
      \ figure0003.jpg

      \ Finished creating file for FRAME = \ 4 \ \ filename =
      \ figure0004.jpg

      \ Finished creating file for FRAME = \ 5 \ \ filename =
      \ figure0005.jpg

      \ Finished creating file for FRAME = \ 6 \ \ filename =
      \ figure0006.jpg

      \ Finished creating file for FRAME = \ 7 \ \ filename =
      \ figure0007.jpg

      \ Finished creating file for FRAME = \ 8 \ \ filename =
      \ figure0008.jpg

      \ Finished creating file for FRAME = \ 9 \ \ filename =
      \ figure0009.jpg

      \ Finished creating file for FRAME = \ 10 \ \ filename =
      \ figure0010.jpg
    </unfolded-io>

    <\input|Python] >
      \;
    </input>

    \;
  </session>

  <section|Sample Output (JPEG)>

  <image|figure0010.jpg||||>
</body>

<\references>
  <\collection>
    <associate|QinitFunc.cpp|<tuple|2|?>>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|3|2>>
    <associate|auto-4|<tuple|4|2>>
    <associate|auto-5|<tuple|5|2>>
    <associate|auto-6|<tuple|6|2>>
    <associate|auto-7|<tuple|7|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Mathematics>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Editing
      Basic Application Files> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Compilation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Editing
      parameter file> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Running
      the example> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Plotting
      in Python> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Sample
      Output (JPEG)> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>
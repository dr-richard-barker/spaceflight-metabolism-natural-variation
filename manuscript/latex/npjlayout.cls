%% npjlayout.cls -- a two-column preprint layout modelled on the 2026 npj Microgravity Article PDF:
%% sans-serif title block and full-width abstract, two-column serif body, unnumbered sans-serif
%% headings, "Fig. 1 |" and "Table 1 |" captions, sans-serif tables, affiliations as a first-page
%% footnote, superscript numbered citations.
%%
%% It deliberately carries NO journal identity: no wordmark, no "Article in Press", no DOI or
%% citation line, no open-access or "Check for updates" badges. A paper set in it is a preprint,
%% and its masthead and running feet say so (\manuscriptstatus). An earlier template in this
%% portfolio stamped a journal citation and DOI on every page of an unpublished manuscript; this
%% class exists so that the layout can be had without that.
%%
%% It provides the sn-jnl commands main.tex uses (\title[short]{}, \author*[n]{}, \fnm, \sur,
%% \email, \orcid, \affil[n]{}, \orgname, \city, \country, \abstract{}, \keywords{}, \backmatter, \bmhead,
%% \botrule), so the source still builds with sn-jnl.cls by changing the \documentclass line.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{npjlayout}[2026/11/05 two-column preprint layout after npj Microgravity Articles]
\DeclareOption*{}% sn-jnl options (pdflatex, sn-nature) are accepted and ignored
\ProcessOptions\relax
\LoadClass[twocolumn,10pt]{article}

%% ---- page, fonts, colour -------------------------------------------------------------------
\RequirePackage[paperwidth=210mm,paperheight=279mm,left=15mm,right=15mm,top=20mm,bottom=20mm,
  headheight=12pt,headsep=7mm,footskip=9mm,columnsep=6mm]{geometry}
\RequirePackage[T1]{fontenc}
\RequirePackage{mathptmx}% Times text and maths; works with amsmath/amssymb loaded afterwards
\RequirePackage[scale=0.92]{tgheros}% Helvetica-like sans for headings, abstract, captions, tables
\RequirePackage{microtype}
\RequirePackage{xcolor}
\definecolor{npjlaccent}{HTML}{2B5C8A}% our own accent; not the journal's brand colour
\definecolor{npjlgrey}{gray}{0.40}
\RequirePackage{marvosym}% \Letter for the corresponding author
\RequirePackage{etoolbox}
\RequirePackage{booktabs}
\RequirePackage{ragged2e}

\renewcommand\normalsize{\@setfontsize\normalsize{9.5}{12}%
  \abovedisplayskip 6pt plus 2pt minus 2pt\belowdisplayskip\abovedisplayskip
  \abovedisplayshortskip 3pt\belowdisplayshortskip 3pt}
\normalsize
\setlength\parindent{1em}
\setlength\parskip{0pt}
\raggedbottom
\setlength\emergencystretch{2em}

%% ---- headings ------------------------------------------------------------------------------
\RequirePackage{titlesec}
\setcounter{secnumdepth}{0}
\titleformat{\section}{\sffamily\bfseries\fontsize{11.5}{13.5}\selectfont}{}{0pt}{}
\titlespacing*{\section}{0pt}{11pt plus 3pt minus 1pt}{4pt}
\titleformat{\subsection}{\sffamily\bfseries\fontsize{9.5}{11.5}\selectfont}{}{0pt}{}
\titlespacing*{\subsection}{0pt}{8pt plus 2pt minus 1pt}{2pt}
\newcommand\backmatter{\par}
\newcommand\bmhead[1]{\par\addvspace{6pt}\noindent{\sffamily\bfseries\fontsize{9.5}{11.5}\selectfont #1}%
  \par\nobreak\addvspace{1.5pt}\noindent\ignorespaces}

%% ---- captions, floats, tables --------------------------------------------------------------
\RequirePackage{caption}
\DeclareCaptionLabelSeparator{npjbar}{\ \textbar\ }
\captionsetup{font={sf,footnotesize}, labelfont=bf, labelsep=npjbar, justification=justified,
  singlelinecheck=false, skip=5pt}
\captionsetup[figure]{name=Fig.}
\captionsetup[table]{name=Table, position=top}
% Floats reset the font when they start, so sans-serif tables are set there, not in a begin hook.
\appto\@floatboxreset{\sffamily\footnotesize}
\setlength\heavyrulewidth{0.6pt}
\setlength\lightrulewidth{0.35pt}
\AtBeginDocument{\providecommand\botrule{\bottomrule}}
\setlength\textfloatsep{12pt plus 2pt minus 2pt}
\setlength\dbltextfloatsep{12pt plus 2pt minus 2pt}
\renewcommand\topfraction{0.9}
\renewcommand\dbltopfraction{0.9}
\renewcommand\textfraction{0.08}
\renewcommand\floatpagefraction{0.8}
\renewcommand\dblfloatpagefraction{0.8}

%% ---- citations ------------------------------------------------------------------------------
\RequirePackage[numbers,super,sort&compress]{natbib}
\bibliographystyle{sn-nature}
\setlength\bibsep{1.5pt}
\renewcommand\bibfont{\footnotesize}
\renewcommand\bibnumfmt[1]{#1.}
\providecommand\bibcommenthead{}
\renewcommand\refname{References}

%% ---- front matter: the sn-jnl commands ------------------------------------------------------
\gdef\@npjshorttitle{}
\gdef\@npjauthors{}
\gdef\@npjaffils{}
\gdef\@npjemail{}
\gdef\@abstract{}
\gdef\@keywords{}
\newif\if@npjfirstauthor \@npjfirstauthortrue
\RenewDocumentCommand\title{o m}{\gdef\@title{#2}\IfValueT{#1}{\gdef\@npjshorttitle{#1}}}
\RenewDocumentCommand\author{s o m}{%
  \if@npjfirstauthor\global\@npjfirstauthorfalse\else\g@addto@macro\@npjauthors{, }\fi
  \g@addto@macro\@npjauthors{#3}%
  \IfValueT{#2}{\g@addto@macro\@npjauthors{\textsuperscript{#2}}}%
  \IfBooleanT{#1}{\g@addto@macro\@npjauthors{\,\textsuperscript{\Letter}}}}
\newcommand\fnm[1]{#1}
\newcommand\sur[1]{#1}
\newcommand\orgname[1]{#1}
\newcommand\city[1]{#1}
\newcommand\country[1]{#1}
\newcommand\email[1]{\gdef\@npjemail{#1}}
\gdef\@npjorcid{}
\newcommand\orcid[1]{\gdef\@npjorcid{#1}}
\NewDocumentCommand\affil{o m}{\g@addto@macro\@npjaffils{\IfValueT{#1}{\textsuperscript{#1}}#2. }}
\renewcommand\abstract[1]{\gdef\@abstract{#1}}
\newcommand\keywords[1]{\gdef\@keywords{#1}}

%% ---- preprint status (shown on every page) --------------------------------------------------
\gdef\@npjstatus{Preprint, not peer reviewed}
\gdef\@npjversion{}
\gdef\@npjurl{}
\newcommand\manuscriptstatus[1]{\gdef\@npjstatus{#1}}
\newcommand\manuscriptversion[1]{\gdef\@npjversion{#1}}
\newcommand\manuscripturl[1]{\gdef\@npjurl{#1}}

%% ---- running heads and feet -----------------------------------------------------------------
\RequirePackage{fancyhdr}
\fancypagestyle{npjlfirst}{\fancyhf{}%
  \fancyfoot[L]{\sffamily\footnotesize\color{npjlaccent}\@npjstatus}%
  \fancyfoot[R]{\sffamily\footnotesize\thepage}%
  \renewcommand\headrulewidth{0pt}\renewcommand\footrulewidth{0pt}}
\fancypagestyle{npjlmain}{\fancyhf{}%
  \fancyhead[L]{\sffamily\footnotesize\color{npjlgrey}\@npjshorttitle}%
  \fancyhead[R]{\sffamily\footnotesize\bfseries Article}%
  \fancyfoot[L]{\sffamily\footnotesize\color{npjlaccent}\@npjstatus}%
  \fancyfoot[R]{\sffamily\footnotesize\thepage}%
  \renewcommand\headrulewidth{0pt}\renewcommand\footrulewidth{0pt}}
\pagestyle{npjlmain}

%% ---- title block ----------------------------------------------------------------------------
\newcommand\@npjtitleblock{%
  \noindent{\sffamily\bfseries\fontsize{15}{17}\selectfont\color{npjlaccent}Preprint}%
  \hfill{\sffamily\bfseries\fontsize{14}{17}\selectfont Article}\par
  \vspace{1pt}%
  \noindent{\sffamily\fontsize{7.5}{9}\selectfont\color{npjlgrey}\@npjstatus
    \ifx\@npjversion\@empty\else\ \textbar\ \@npjversion\fi}%
  \hfill{\sffamily\fontsize{7.5}{9}\selectfont\color{npjlaccent}\@npjurl}\par
  \vspace{16pt}%
  {\RaggedRight\sffamily\bfseries\fontsize{22}{26}\selectfont\@title\par}%
  \vspace{12pt}%
  {\RaggedRight\sffamily\bfseries\fontsize{10}{13}\selectfont\@npjauthors\par}%
  \vspace{10pt}{\color{npjlgrey}\hrule height 0.4pt}\vspace{8pt}%
  {\sffamily\fontsize{10}{13.6}\selectfont\@abstract\par}%
  \ifx\@keywords\@empty\else
    \vspace{6pt}{\sffamily\fontsize{7.5}{9.5}\selectfont\color{npjlgrey}\textbf{Keywords:}
      \@keywords\par}\fi
  \vspace{8pt}{\color{npjlgrey}\hrule height 0.4pt}\vspace{14pt}}
\renewcommand\maketitle{%
  \twocolumn[\@npjtitleblock]%
  \thispagestyle{npjlfirst}%
  {\let\thefootnote\relax
   \footnotetext{\sffamily\fontsize{7}{8.5}\selectfont\@npjaffils
     \ifx\@npjemail\@empty\else\Letter\ e-mail: \@npjemail\fi
     \ifx\@npjorcid\@empty\else; ORCID \@npjorcid\fi}}%
  \setcounter{footnote}{0}}

%% ---- links (after the document's own packages) ----------------------------------------------
\AtEndPreamble{%
  \RequirePackage{hyperref}\urlstyle{same}%
  \hypersetup{colorlinks=true,linkcolor=black,citecolor=npjlaccent,urlcolor=npjlaccent,
    pdfcreator={npjlayout.cls (preprint layout; no journal identity)}}}
\endinput

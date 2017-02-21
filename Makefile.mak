# file:       Makefile.mak
# author:     Andrea Vedaldi (edited for vlg by Taehee Lee)
# descrption: Microsoft NMake makefile

# Customization:
# - MATLABROOT : must point to MATLAB root directory (undef = no MEX support)

# MATLABROOT = C:\MATLAB7
MATLABROOT = "C:\Program Files\MATLAB\R2007b"

# --------------------------------------------------------------------
#                                                                Flags
# --------------------------------------------------------------------
# Debug info is embedded in .obj and .lib files (CodeView /Z7 option)
# but in .pdb files for .exe and .dll (since the linker does not
# produce CodeView output anymore).
#
# CFLAGS
#   /nologo            : CL does not display splash
#   _CRT_NO_DEPRECATE  : Do not deprecate `insecure' fscanf, snprintf, ...
#   __LITTLE_ENDIAN__  : Signal little endian architecture
#   /I.                : Add VLROOT to include search path
#   /Z7                : Embedded CodeView debug info in .obj
#   /MT                : Multi-thread run-time library
#   /TC                : Source code is C (not C++)
#   /W3                : Usa all warnings
#   /Zp8               : Align structures to 8 bytes
#   /O2                : Optimize for speed
#   /arch:SSE2         : Enable SSE2 instructions
#
# LFLAGS
#   /NOLOGO            : LINK does not display splash
#   /INCREMENTAL:NO    : No incremental linking
#   /MANIFEST:NO       : No manifest
#   /DEBUG             : Generate debug info (.pdb files)
#
# MEX_RC               : MEX .rc file location
#
# MEX_CFLAGS
#   /D_WINDLL          : Signal DLL code
#   /DMATLAB_MEX_FILE  : Signal MATLAB MEX code
#
# MEX_LFLAGS
#  /DLL                : Produce a DLL
#  /EXPORT:mexFunction : Export MEX file entry point

vlfeatdir  = ..\vlfeat
bindir     = $(vlfeatdir)\bin\win32

CFLAGS     = /nologo /TC /MT \
             /D"__VISUALC__" /D"WIN32" \
             /D"__LITTLE_ENDIAN__" \
             /D"_CRT_SECURE_NO_DEPRECATE" \
             /I. /I$(vlfeatdir) \
	     /Itoolbox /I$(vlfeatdir)\toolbox \
             /W1 /Wp64 /Z7 /Zp8 /O2 /arch:SSE2

LFLAGS     = /NOLOGO /INCREMENTAL:NO /MANIFEST:NO \
             /LIBPATH:$(bindir) vl.lib \
             /DEBUG

MEX_RC     = $(MATLABROOT)\extern\include\mexversion.rc

MEX_CFLAGS = $(CFLAGS) /I$(MATLABROOT)\extern\include \
             /DMATLAB_MEX_FILE /D_WINDLL

MEX_LFLAGS = $(LFLAGS) \
             /DLL /EXPORT:mexFunction \
             /MACHINE:X86 \
             /LIBPATH:$(MATLABROOT)\extern\lib\win32\microsoft \
             libmx.lib libmex.lib libmat.lib  


mexdir = \
 toolbox\bundle

mexsrc =                                   \
 toolbox\bundle\mex_bundle_1_XABeUVWeAeB.c \
 toolbox\bundle\mex_bundle_2_Se_.c         \
 toolbox\bundle\mex_bundle_3_db_new.c      \
 toolbox\bundle\mex_bundle_proj_1_XABeUVWeAeB.c \
 toolbox\bundle\mex_bundle_proj_2_Se_.c         \
 toolbox\bundle\mex_bundle_proj_3_db_new.c

mexdll = $(mexsrc:.c=.dll)
mexres = $(mexsrc:.c=.res)
mexpdb = $(mexsrc:.c=.pdb)

all: $(mexdll)

# --------------------------------------------------------------------
#                                                    Maintenance rules
# --------------------------------------------------------------------

clean:
	-del $(mexpdb)

distclean: clean
	-del $(mexdll)
info:
	@echo ** mexsrc     = $(mexsrc)
	@echo ** mexdll     = $(mexdll)
	@echo ** CC         = $(CC)
	@echo ** CFLAGS     = $(CFLAGS)
	@echo ** MEX_CFLAGS = $(MEX_CFLAGS)
	@echo ** MEX_LFLAGS = $(MEX_LFLAGS)

# --------------------------------------------------------------------
#                                                          Build rules
# --------------------------------------------------------------------

# *.c -> *.dll
{$(mexdir)}.c{$(mexdir)}.dll:
	@echo CC  $(<) ===^> $(@R).dll
	@$(CC) $(MEX_CFLAGS) /c /Fo"$(@R).obj" "$(<)"
	@RC /fo"$(@R).res" $(MEX_RC)
	@LINK $(MEX_LFLAGS) "$(@R).res" "$(@R).obj" /OUT:$(@)
	@-del "$(@R).obj"
	@-del "$(@R).exp"
	@-del "$(@R).lib"
	@-del "$(@R).res"


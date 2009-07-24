function _c3d() 
{
  local cur prev opts interp
  opts=""
  opts="$opts -add"
  opts="$opts -anisotropic-diffusion  -ad"
  opts="$opts -antialias  -alias"
  opts="$opts -as  -set"
  opts="$opts -background"
  opts="$opts -binarize"
  opts="$opts -centroid"  
  opts="$opts -connected-components  -connected  -comp"
  opts="$opts -clear"
  opts="$opts -clip"
  opts="$opts -copy-transform  -ct"
  opts="$opts -create"
  opts="$opts -dilate"
  opts="$opts -divide"
  opts="$opts -dup -duplicate"
  opts="$opts -endfor"
  opts="$opts -erode"
  opts="$opts -erf"
  opts="$opts -exp"
  opts="$opts -fft"
  opts="$opts -foreach"
  opts="$opts -glm"
  opts="$opts -histmatch  -histogram-match"
  opts="$opts -info"
  opts="$opts -info-full"
  opts="$opts -insert  -ins"
  opts="$opts -interpolation  -interp  -int"
  opts="$opts -iterations"
  opts="$opts -label-statistics  -lstat"
  opts="$opts -laplacian  -laplace"
  opts="$opts -levelset"
  opts="$opts -levelset-curvature"
  opts="$opts -levelset-advection"
  opts="$opts -ln  -log"
  opts="$opts -log10"
  opts="$opts -mcs -multicomponent-split"
  opts="$opts -mean"
  opts="$opts -merge"
  opts="$opts -mi  -mutual-info"
  opts="$opts -mixture  -mixture-model"
  opts="$opts -multiply  -times"
  opts="$opts -nmi  -normalized-mutual-info"
  opts="$opts -normpdf"
  opts="$opts -noround"
  opts="$opts -nospm"
  opts="$opts -o"
  opts="$opts -omc -output-multicomponent"
  opts="$opts -oo -output-multiple"
  opts="$opts -origin"
  opts="$opts -overlap"
  opts="$opts -pad"
  opts="$opts -pixel"
  opts="$opts -pop"
  opts="$opts -popas"
  opts="$opts -probe"
  opts="$opts -push  -get"
  opts="$opts -reciprocal"
  opts="$opts -region"
  opts="$opts -replace"
  opts="$opts -resample"
  opts="$opts -resample-mm"
  opts="$opts -reslice-itk"
  opts="$opts -reslice-matrix"
  opts="$opts -reslice-identity"
  opts="$opts -rms"
  opts="$opts -round"
  opts="$opts -scale"
  opts="$opts -shift"
  opts="$opts -signed-distance-transform  -sdt"
  opts="$opts -smooth"
  opts="$opts -spacing"
  opts="$opts -split"
  opts="$opts -sqrt"
  opts="$opts -staple"
  opts="$opts -spm"
  opts="$opts -stretch"
  opts="$opts -threshold  -thresh"
  opts="$opts -trim"
  opts="$opts -trim-to-size"
  opts="$opts -type"
  opts="$opts -verbose"
  opts="$opts -version"
  opts="$opts -vote"
  opts="$opts -vote-label"
  opts="$opts -voxel-sum"
  opts="$opts -voxel-integral  -voxel-int"
  opts="$opts -warp"
  opts="$opts -warp-label  -warplabel  -wl"

  interp="nearest linear cubic sinc gaussian Nearest Linear Cubic Sinc Gaussian"

  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  if [[ ${cur} == -* ]] ; then
    COMPREPLY=( `compgen -W "${opts}" -- ${cur}` )
    return 0

  elif [[ ${prev} == -int* ]] ; then
    COMPREPLY=( `compgen -W "${interp}" -- ${cur}` )
    return 0

  else
    COMPREPLY=( `compgen -f ${cur}` )
    return 0
  fi
}
complete -F _c3d -o filenames c3d


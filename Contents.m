% HYBRIDMT
%
% Files
%   drawhudsonnet     - Draw Hudson network and/or provide MT projection to [u,v] space.
%   drawstereonet     - Create stereonet plot.
%   focimt            - Perform seismic full moment tensor inversion using fociMT software.
%   genmt_raw         - GETMT_RAW Generate synthetic event/phase data for fociMT in RAW ASCII format.
%   genmt_vel1d       - Generate synthetic event/phase data for fociMT in 1D velocity model ASCII format
%   getsolution       - Extract parameter(s) from moment tensor solution cell array
%   hybridmt          - Perform refinement of moment tensors using hybridMT technique.
%   hybridmt_install  - Add current directory to MATLAB path and build searcheable database.
%   hybridmt_show     - Internal function to generate figures from hybridMT refinement.
%   project_stereonet - Project data onto stereonet plot.
%   readraw           - Read fociMT input ASCII file in RAW format.
%   readsolution      - Read seismic moment tensor solution file created by fociMT
%   readvel1d         - Read fociMT input ASCII file in 1D velocity model format.
%   rpgen             - Calculate radiation pattern using shear-tensile source model.
%   writeinput        - Converts ASCII input file in RAW or 1D velocity format to input cell array

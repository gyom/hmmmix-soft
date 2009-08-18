% Summer 2009. Guillaume Alain, gyomalin@gmail.com


%%% Installation procedure. %%%


% For the sake of the discussion, let's assume that you want to install the
% code for hmmmix-soft code into the directory
% /home/gyom/.matlab/work/hmmsoft
% and that the contents of the zip archive has alread been expanded into that
% directory, keeping the internal structure of the archive.
% I don't personally recommend adding anything other than the main directory into the
% Matlab path. That is, don't add it recursively because you will catch the 'optional_kmedoids_wrapper'
% directory that overshadows some of Sohrab Shah's functions from his CNAhmmer package.



cd /home/gyom/.matlab/work/hmmsoft

% Compile the C sources for certain functions. If you get a problem, open
% the compile_hmmmix_soft.m file and run the commands line by line to diagnose the problem.
compile_hmmmix_soft

% Now let's run some unit tests to be sure that some of the basic functions work properly.
% In the 'reference' directory there are Matlab implementations for functions such as
% the forwards-backwards algorithm whose results should match those of the C implementation.
% You can always add the 'reference' directory to your path if you don't mind catching these
% functions, but they are not required for the rest of the algorithm. I would advise to only
% cd into that directory at installation to check the C functions compiled and then never go
% into that directory again.
cd reference
unit_test_normalize
unit_test_viterbi_path
unit_test_fwd_back_MatlabC
unit_test_hmmmix_frugal_hM_KTg_MatlabC



% Inside the 'reference' directory, you can run this function to make sure that both
% implementations of hmmmix-soft (in memory and on hard drive) produce the same results.
% This is more of a sanity check than anything. If this passes, it probably means that
% you're clear and there is not problem in the installation. This script can also serve as
% an example of usage for my implementation of the hmmmix-soft algorithm.

script_to_compare_hmmsoft_inRAM_vs_onHD

% The only potential problem now is if you want to use my function that calls
% Sohrab's k-medoids code. If you want to do so, you need to add the 'optional_kmedoids_wrapper'
% directory to the path. I don't really like recommending using this because I'm polluting the
% namespace by overshadowing some of Sohrab's functions, but in some cases it might be what you
% want to do. If you're using Sohrab's CNAmmer code too, you can always add the
% 'optional_kmedoids_wrapper' directory only when you need it without saving the updated Matlab path.

 
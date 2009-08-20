%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
% UBC 2007-2009. Guillaume Alain, gyomalin@gmail.com    %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Some notes
% 2. Installation procedure.
% 3. About k-medoids
% 4. Permission to use the code.

%%%%%%%%%%%%%%%%%%%%%%
%%% 1. Some notes  %%%
%%%%%%%%%%%%%%%%%%%%%%

% This code contains the essential pieces of the hmmmix-soft algorithm described in my
% thesis at UBC, entitled "Model-based clustering for aCGH data using variational EM".
% All the scripts to orchestrate the experiences that I carried are not found here, but
% there is one file called "script_to_compare_hmmsoft_inRAM_vs_onHD.m" that should serve
% as a basis for using this code.
%
% The hmmmix-soft algorithm is my version of the hmmmix algorithm from Sohrab Shah. Refer
% to Sohrab Shah's PhD thesis if needed.
%
% This code is provided more or less "as is", but I did try to do a good job at packaging
% it so I could be built from scratch and used easily. My apologies to anyone finding the
% documentation insufficient. The final steps in the production of this "package" for
% distribution to the general public were done after I was finished with my experiments
% for my thesis. This means that, although I have tested my code, performed unit tests and
% such, I've never really used the final package myself. This limits how confident I can be
% when I say that it's working fine.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2. Installation procedure  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This procedure was tested on Mac OSX 10.5.8 with Matlab R2008b
% and on Linux (Open Suse 10.3) with Matlab R2007b.
% To install this on Windows you need a C compiler because a lot of functions
% were written in C when vectorization was impossible. There is nothing
% fundamentally related to Linux/Unix in the code that would cause problems,
% but the compilation method might differ slightly depending on your compiler.
%
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
% an example of usage for my implementation of the hmmmix-soft algorithm. There is a small issue
% where the results can differ by almost nothing. I decided to draw the line there and not
% track down this minor issue. See inside the script for further comments.

script_to_compare_hmmsoft_inRAM_vs_onHD
cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 3. About k-medoids  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There is a function called 'kmedoids_inRAM_wrapper_to_Sohrab_code' that uses
% Sohrab Shah's k-medoids code to initialize the patients assignments. To use it, you
% need to have his CNA-HMMer code on your Matlab path. On top of that, you need to add the 
% 'hmmmix_soft/optional_kmedoids_wrapper' directory to your path when using the function
% 'kmedoids_inRAM_wrapper_to_Sohrab_code'. It was necessary for me to overshadow some of Sohrab's
% original functions to get the functionality that I wanted, but this is not something that you
% want to do permanently if you are using his CNA-HMMer code regularly.
%
% I'm personally using a more recent version of the code available on his web site. The code that I
% use dates from around January 2009 and the directory is named "CNA-HMMer-spec". I just added 
% all the directory from his code recursively and it works for me, but the important thing is that
% the functions from my 'optional_kmedoids_wrapper' are found BEFORE those of CNA-HMMer-spec that they
% are intended to overshadow.
%
% You can test the k-medoids functions by uncommenting lines 67-68 in
% 'script_to_compare_hmmsoft_inRAM_vs_onHD'. This should indicate how the 
% 'kmedoids_inRAM_wrapper_to_Sohrab_code' function is meant to be used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4. Permission to use the code  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I don't think there is much money to be made with this code, but I'd still like to get credit
% for what I did. So, there you have it. That makes it more of a BSD license than a GPL license.
% I would be really happy if my work was actually used by others.

% I wrote everything in the root directory except mexutils.c, mexutils.h, repmatC.c and normalize.m.
% Some functions in the 'reference' directory come from Kevin Murphy and Sohrab Shah. Most of what is
% in the 'option_kmedoids_wrapper' comes from Sohrab Shah more or less directly.

% Finally, I should say that, at the time I release this code, I don't plan to work on it anymore.
% It was a nice thing to be able to work on this problem for my thesis, but it's not really my "field"
% and I won't be the one maintaining this code if it's actually used for serious research. In fact,
% if you start finding bugs that you can fix and if you want to integrate this code into your own,
% you are free to do so. It'd be nice to send me an email if you do that, though.



 
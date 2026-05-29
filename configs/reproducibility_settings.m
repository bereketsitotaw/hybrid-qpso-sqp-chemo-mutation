function repro = reproducibility_settings()
%REPRODUCIBILITY_SETTINGS Pinned RNG and solver options for reproduction.

    repro = struct();
    repro.matlab_version = 'R2024a';
    repro.rng_seed = 42;
    repro.rng_generator = 'twister';
    repro.N = 100;
end

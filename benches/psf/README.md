# PSF Tests

This subrepository is made to be able to evaluate the performance of the PSF.
For large parameter sets, the time to generate the keys increases dramatically.
To only evaluate the [`sample_p`] algorithm, we will use the code of this library to generate an instantation of the scheme that can be loaded in the dedicated tests without having to generate the keys itself.

"""
    plot_input_pdfs(
        pdf_params::Union{BernsteinPDFParams,BernsteinDirichletPDFParams};
        xmin::Real, xmax::Real, nx::Integer
    )

Plot PDFs using the given PDF parameters.

!!! Note
    Package `Plots` must be loaded, e.g. via `import Plots`, to use this function.
"""
function plot_input_pdfs end
export plot_input_pdfs



"""
    plot_model_space(pdf_params, samples)

Compare truth and posterior samples in the model space.

!!! Note
    Package `Plots` must be loaded, e.g. via `import Plots`, to use this function.
"""
function plot_model_space end
export plot_model_space


"""
    plot_data_space(
        pdf_params, sim_data, samples,
        qcdnum_params, splint_params, quark_coeffs, md
    )

Compare truth and posterior samples in the data space.

!!! Note
    Package `Plots` must be loaded, e.g. via `import Plots`, to use this function.
"""
function plot_data_space end
export plot_data_space

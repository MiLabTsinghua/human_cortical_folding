def do_perturb(oracle,TF,n_grid = 40,plot = False):
    oracle.simulate_shift(perturb_condition={TF: 0.0},
                          n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=200,
                                    knn_random=True,
                                    sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    if plot:
        fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
    
        scale = 50
        Show quiver plot
        oracle.plot_quiver(scale=scale, ax=ax[0])
        ax[0].set_title(f"Simulated cell identity shift vector: {TF} KO")
        
        Show quiver plot that was calculated with randomized graph.
        oracle.plot_quiver_random(scale=scale, ax=ax[1])
        ax[1].set_title(f"Randomized simulation vector")
        
        plt.show()

    
    oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    print(TF,flush=True)
    oracle.suggest_mass_thresholds(n_suggestion=12)
    print('*'*50,flush=True)
    return(oracle)



    
def calculate_PS(oracle,gradient,TF):
    from celloracle.applications import Oracle_development_module
    dev = Oracle_development_module()

    # Load development flow
    dev.load_differentiation_reference_data(gradient_object=gradient)
    
    # Load simulation result
    dev.load_perturb_simulation_data(oracle_object=oracle)
    
    
    # Calculate inner produc scores
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)

    dev.inner_product_df.to_csv(f'Perturbation_Score/{TF}.csv')
    return(dev)






def prepare_gradient(oracle,pseudotime_key,min_mass = 18,n_grid = 40):
    from celloracle.applications import Gradient_calculator
    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key=pseudotime_key)
    gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    gradient.calculate_mass_filter(min_mass=min_mass, plot=False)
    gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=False)
    gradient.calculate_gradient()
    return(gradient)



def draw_TF_perturabation(oracle,TF,min_mass,save=True):
    oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
    fig, ax = plt.subplots(figsize=[6, 6])
    oracle.plot_cluster_whole(ax=ax, s=10)
    scale_simulation = 8
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
    ax_cbar = fig.add_axes([0.2, 0.3, 0.03, 0.6])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # dummy array required for ScalarMappable
    
    # 添加 colorbar legend
    cbar = plt.colorbar(sm, cax=ax_cbar)
    cbar.set_label("Pseudoage Score")
    
    ax.set_title(f'{TF}_perturbation')
    if save == True:
        plt.savefig(f'TF_perturbation/output_celltype_sample/figures/{TF}.pdf')
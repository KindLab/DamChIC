def simulate_dynamics(n_obs, switch, direction=1, noise_level=.1, seed=0):
    
    """
    Based on: `scvelo.datasets.simulation` 10.1038/s41587-020-0591-3 
     
    Simulation of dynamics at specified switch. 
    
    The parameters for each reaction are randomly sampled from a log-normal distribution
    and time events follow the Poisson law.
    
    Parameters
    ----------
    n_obs : int, number of observations.
    
    switch : float between 0 and 1. Fraction of n_obs at which dynamic starts.
    
    direction : int 1 or -1. 1: upregulated, -1: downregulated dynamic. 
    
    noise_level : float between 0 and 1. 
    
    """
    np.random.seed(seed)
    switch = int(switch * n_obs)
    
    def draw_poisson(n_obs):
        t = np.cumsum([-0.1 * np.log(np.random.uniform(0, 1)) for _ in range(n_obs - 1)])
        return np.insert(t, 0, 0) # prepend t0=0
    
    t = draw_poisson(n_obs)
    t *= direction
    t_switch = np.insert(t, 0, np.repeat(0,switch))[:-switch]
    t_switch = (t_switch - t_switch.min()) / (t_switch.max() - t_switch.min())
    
    def draw_normal(t_switch, noise_level=noise_level):
        return np.hstack([np.random.normal(v, scale=noise_level) for v in t_switch])

    t_normal = draw_normal(t_switch)
    t_normal = (t_normal - t_normal.min()) / (t_normal.max() - t_normal.min())
    
    return t_normal

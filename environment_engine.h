#ifndef ENVIRONMENT_ENGINE_H
#define ENVIRONMENT_ENGINE_H

#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <map>
#include <any>
#include <optional>
#include <functional>
#include <memory>


namespace EnvEngine {

    // Precision type as specified
    using real = double;

    // Data type for a 3D vector
    struct Vector3 {
        real x = 0.0, y = 0.0, z = 0.0;
    };

    // Global physical constants
    namespace GlobalConstants {
        constexpr real G = 9.80665;
        constexpr real R_d = 287.058;
        constexpr real R_v = 461.495;
        constexpr real GAMMA_AIR = 1.4;
        constexpr real SIGMA_SB = 5.670374419e-8;
        constexpr real LATENT_HEAT_VAPORIZATION_LV = 2.501e6;
        constexpr real SPECIFIC_HEAT_AIR_CP = 1005.0;
        constexpr real WATER_DENSITY = 1000.0;
        constexpr real PI = 3.14159265358979323846;
        constexpr real SOLAR_CONSTANT = 1361; // W/m^2
        constexpr real VON_KARMAN_K = 0.41;
        // FINAL FIX: Corrected evaporation coefficient based on physical constants
        constexpr real EVAPORATION_COEFF = 2.5e-7; 
    }
    
    // --- Utility Classes ---
    class PerlinNoise {
    public:
        PerlinNoise(unsigned int seed);
        real noise(real x, real y, real z) const;
    private:
        std::vector<int> p;
        real fade(real t) const;
        real lerp(real t, real a, real b) const;
        real grad(int hash, real x, real y, real z) const;
    };

    // --- Forward Declarations ---
    class WindAndTurbulence;
    class SurfaceThermodynamics;
    class SolarRadiation;
    class AtmosphereCore;
    class PrecipitationAndHydrometeors;
    class MoistureAndHydrology;
    class Environment;

    // --- Module Declarations ---

    class TimeManager {
    public:
        void advance(real dt_s);
        void set_time(const std::string& iso_string);
        Vector3 compute_sun_vector(real latitude_deg, real longitude_deg) const;
        real get_sim_time() const { return sim_time_s; }
        std::string get_utc_time_str() const;
    private:
        std::chrono::system_clock::time_point utc_time;
        std::tm current_utc_tm;
        real sim_time_s = 0.0;
    };

    class AtmosphereCore {
    public:
        struct State {
            real T_air_C, p_air_kPa, RH_percent, rho_air_kg_m3;
        };
        State query_state_at(const Vector3& position) const;
    private:
        real sea_level_T_air_C = 25.0, sea_level_p_air_kPa = 101.325, sea_level_RH_percent = 75.0, lapse_rate_C_per_m = -0.0065;
    };
    
    class StochasticEventManager {
    public:
        struct Event {
            std::string type;
            std::map<std::string, real> params;
            real start_time, duration;
            Vector3 position, velocity;
        };
        void initialize(unsigned int seed);
        void update(real dt_s, real current_sim_time, const WindAndTurbulence& wind_module);
        void add_event(const Event& new_event) { active_events.push_back(new_event); }
        const std::vector<Event>& get_active_events() const { return active_events; }
    private:
        std::mt19937 rng;
        std::vector<Event> active_events;
        std::uniform_real_distribution<real> uniform_dist{0.0, 1.0};
        real time_since_last_event = 0.0;
        real poisson_rate_events_per_second = 1.0 / 1800.0;
        std::discrete_distribution<> event_type_dist;
    };
    
    class CloudModel {
    public:
        void initialize(unsigned int seed);
        void update(real dt_s, const Vector3& wind_advection);
        real get_cloud_cover_at(const Vector3& position, real sim_time) const;
    private:
        std::unique_ptr<PerlinNoise> cloud_noise;
        Vector3 advection_offset = {0,0,0};
    };

    struct TrackGeometry {
        std::vector<Vector3> shadow_casters_min;
        std::vector<Vector3> shadow_casters_max;
        std::vector<real>    caster_heights;
    };

    class SolarRadiation {
    public:
        struct Irradiance { real I_direct_W_m2, I_diffuse_W_m2, I_longwave_down_W_m2, zenith_angle_deg; bool is_in_shadow; };
        Irradiance query_irradiance(const Vector3& position, const Vector3& sun_vector, const AtmosphereCore::State& air_state, real local_cloud_cover, const TrackGeometry& track) const;
    private:
        real aerosol_optical_depth = 0.2;
    };

    class WindAndTurbulence {
    public:
        struct WindState { Vector3 wind_vector_m_s; real turbulence_sigma_u; };
        void initialize(unsigned int seed);
        void set_base_wind(real speed_ms, real direction_deg) {
            ref_wind_speed.x = speed_ms * cos(direction_deg * GlobalConstants::PI / 180.0);
            ref_wind_speed.y = speed_ms * sin(direction_deg * GlobalConstants::PI / 180.0);
        }
        WindState query_wind_at(const Vector3& position, real sim_time_s, const std::vector<StochasticEventManager::Event>& events) const;
        Vector3 get_mean_wind_at_ref_height() const { return ref_wind_speed; };
    private:
        Vector3 ref_wind_speed = {4.0, 1.5, 0.0};
        real ref_height = 10.0, surface_roughness_z0 = 0.03, friction_velocity_u_star = 0.3;
        std::unique_ptr<PerlinNoise> turbulence_noise;
    };
    
    class PrecipitationAndHydrometeors {
    public:
        struct PrecipData { 
            real rain_intensity_mm_hr = 0.0;
            real mean_drop_diameter_mm = 0.0;
            real dsd_shape_param = 0.0;
        };
        PrecipData query_rain_at(const Vector3& position, const std::vector<StochasticEventManager::Event>& events) const;
    };

    class MoistureAndHydrology {
    public:
        struct MoistureFlux { real precipitation_mm_s, evaporation_mm_s, runoff_mm_s; };
        void update(real dt_s, const PrecipitationAndHydrometeors::PrecipData& precip, const AtmosphereCore::State& air, 
                    const SurfaceThermodynamics& surface, const WindAndTurbulence::WindState& wind, real terrain_slope);
        
        void add_water_flow(real flow_mm) { surface_water_depth_mm += flow_mm; }
        MoistureFlux get_fluxes() const { return current_flux; }
        real get_surface_water_depth_mm() const { return surface_water_depth_mm; }
    private:
        MoistureFlux current_flux = {0.0, 0.0, 0.0};
        real surface_water_depth_mm = 0.0;
        real soil_moisture_content = 0.5;
    };

    class SurfaceThermodynamics {
    public:
        void update(real dt_s, const SolarRadiation::Irradiance& solar, const AtmosphereCore::State& air, const MoistureAndHydrology& hydro);
        real get_surface_temp_C() const { return surface_temp_C; }
    private:
        real surface_temp_C = 28.0;
        real albedo = 0.12, emissivity = 0.95, heat_capacity_J_m2_K = 1.5e6;
    };

    class AerosolsAndVisibility {
    public:
        void update(real dt_s, const WindAndTurbulence::WindState& wind, const PrecipitationAndHydrometeors::PrecipData& precip);
        real get_visibility_m(const Vector3& position, real rh_percent, const std::vector<StochasticEventManager::Event>& events) const;
    private:
        real aerosol_conc_ug_m3 = 40.0;
    };
    
    class SurfaceGrid {
    public:
        // FIX: Moved Cell struct to public so main.cpp can access its members
        struct Cell {
            MoistureAndHydrology hydro_module;
            SurfaceThermodynamics thermo_module;
            Vector3 center_pos;
            real elevation_m = 0.0;
            real terrain_slope = 0.01;
        };

        void initialize(int x_dim, int y_dim, real cell_size_m, const Vector3& origin);
        void update_all(real dt_s, const Environment& env);
        
        const Cell& get_cell_at(const Vector3& position) const;

        int get_x_dim() const { return grid_x_dim; }
        int get_y_dim() const { return grid_y_dim; }
        real get_cell_size() const { return cell_size; }
        // FIX: Added a public accessor for the grid origin
        const Vector3& get_origin() const { return grid_origin; }
        
        int get_x_idx(real world_x) const;
        int get_y_idx(real world_y) const;

    private:
        std::vector<std::vector<Cell>> grid;
        int grid_x_dim = 0, grid_y_dim = 0;
        real cell_size = 0.0;
        Vector3 grid_origin;
    };

    struct VehicleState {
        Vector3 position;
        real speed_ms;
    };

    class Environment {
    public:
        // --- Decoupled Track State ---
        struct TrackConditionRate {
            real rubbering_in_rate; // kg/m^2/s
            real marble_generation_rate; // kg/m^2/s
        };

        struct AtmosphereQueryResult { real T_air_C, p_air_kPa, RH_percent; Vector3 wind_vector_m_s; real visibility_m; };
        struct SolarQueryResult { real I_direct_W_m2, I_diffuse_W_m2, zenith_angle_deg; bool is_in_shadow; };
        using PrecipitationQueryResult = PrecipitationAndHydrometeors::PrecipData;
        struct SurfaceStateResult { real ground_temp_C, water_depth_mm, rubber_level; }; // No longer managed here
        struct SurfaceModifierResult { 
            real base_grip_modifier;
            real rolling_resistance_delta;
        };
        
        void initialize(const std::map<std::string, std::any>& config, TrackGeometry track);
        void step(real dt_s, const std::vector<VehicleState>& vehicles);

        AtmosphereQueryResult queryAtmosphere(const Vector3& position) const;
        SolarQueryResult querySolar(const Vector3& position) const;
        PrecipitationQueryResult queryPrecipitation(const Vector3& position) const;
        SurfaceStateResult querySurfaceState(const Vector3& position) const;
        SurfaceModifierResult querySurfaceModifiers(const Vector3& position, real vehicle_speed_ms) const;
        const std::vector<std::vector<TrackConditionRate>>& get_track_condition_rates() const;

        const TimeManager& get_time_manager() const { return time_manager; }
        const AtmosphereCore& get_atmosphere_core() const { return atmosphere_core; }
        const SolarRadiation& get_solar_radiation() const { return solar_radiation; }
        const WindAndTurbulence& get_wind_module() const { return wind_and_turbulence; }
        const StochasticEventManager& get_event_manager() const { return stochastic_event_manager; }
        StochasticEventManager& get_event_manager_non_const() { return stochastic_event_manager; }
        const PrecipitationAndHydrometeors& get_precipitation_module() const { return precipitation_module; }
        const CloudModel& get_cloud_model() const { return cloud_model; }
        const SurfaceGrid& get_surface_grid() const { return surface_grid; }
        real get_latitude() const { return latitude_deg; }
        real get_longitude() const { return longitude_deg; }

    private:
        TimeManager time_manager;
        AtmosphereCore atmosphere_core;
        SolarRadiation solar_radiation;
        WindAndTurbulence wind_and_turbulence;
        StochasticEventManager stochastic_event_manager;
        PrecipitationAndHydrometeors precipitation_module;
        AerosolsAndVisibility visibility_module;
        CloudModel cloud_model;
        SurfaceGrid surface_grid;
        TrackGeometry track_geometry;

        std::vector<std::vector<TrackConditionRate>> track_condition_rates;
        
        bool initialized = false;
        real latitude_deg = 12.01, longitude_deg = 79.86;
    };
}
#endif


#include "environment_engine.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <sstream>
#include <numeric>

#if defined(_WIN32) || defined(_WIN64)
#include <time.h>

time_t timegm(struct tm *tm) { return _mkgmtime(tm); }
#endif

namespace EnvEngine {
    // --- PerlinNoise Implementation ---
    PerlinNoise::PerlinNoise(unsigned int seed) {
        p.resize(256);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), std::mt19937(seed));
        p.insert(p.end(), p.begin(), p.end());
    }
    real PerlinNoise::fade(real t) const { return t * t * t * (t * (t * 6 - 15) + 10); }
    real PerlinNoise::lerp(real t, real a, real b) const { return a + t * (b - a); }
    real PerlinNoise::grad(int hash, real x, real y, real z) const {
        int h = hash & 15;
        real u = h < 8 ? x : y;
        real v = h < 4 ? y : h == 12 || h == 14 ? x : z;
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }
    real PerlinNoise::noise(real x, real y, real z) const {
        int X = (int)floor(x) & 255, Y = (int)floor(y) & 255, Z = (int)floor(z) & 255;
        x -= floor(x); y -= floor(y); z -= floor(z);
        real u = fade(x), v = fade(y), w = fade(z);
        int A = p[X] + Y, AA = p[A] + Z, AB = p[A + 1] + Z;
        int B = p[X + 1] + Y, BA = p[B] + Z, BB = p[B + 1] + Z;
        return lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z), grad(p[BA], x-1, y, z)), lerp(u, grad(p[AB], x, y-1, z), grad(p[BB], x-1, y-1, z))),
                       lerp(v, lerp(u, grad(p[AA+1], x, y, z-1), grad(p[BA+1], x-1, y, z-1)), lerp(u, grad(p[AB+1], x, y-1, z-1), grad(p[BB+1], x-1, y-1, z-1))));
    }

    // --- TimeManager Implementation ---
    void TimeManager::advance(real dt_s) {
        sim_time_s += dt_s;
        utc_time += std::chrono::milliseconds(static_cast<long long>(dt_s * 1000.0));
        time_t tt = std::chrono::system_clock::to_time_t(utc_time);
        current_utc_tm = *gmtime(&tt);
    }
    void TimeManager::set_time(const std::string& iso_string) {
        std::tm tm = {};
        std::stringstream ss(iso_string);
        ss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
        if (ss.fail()) {
             throw std::runtime_error("Failed to parse ISO time string: " + iso_string);
        }
        time_t time = timegm(&tm);
        utc_time = std::chrono::system_clock::from_time_t(time);
        sim_time_s = 0.0;
        time_t tt = std::chrono::system_clock::to_time_t(utc_time);
        current_utc_tm = *gmtime(&tt);
    }
    std::string TimeManager::get_utc_time_str() const {
        std::stringstream ss;
        ss << std::put_time(&current_utc_tm, "%Y-%m-%d %H:%M:%S Z");
        return ss.str();
    }
    Vector3 TimeManager::compute_sun_vector(real latitude_deg, real longitude_deg) const {
        real day_of_year = current_utc_tm.tm_yday + 1;
        real utc_hour = current_utc_tm.tm_hour + current_utc_tm.tm_min / 60.0 + current_utc_tm.tm_sec / 3600.0;
        real phi_rad = latitude_deg * GlobalConstants::PI / 180.0;
        real B_rad = (360.0 / 365.0) * (day_of_year - 81) * GlobalConstants::PI / 180.0;
        real eot_minutes = 9.87 * sin(2 * B_rad) - 7.53 * cos(B_rad) - 1.5 * sin(B_rad);
        real local_solar_time = utc_hour + (longitude_deg / 15.0) + (eot_minutes / 60.0);
        real hour_angle_H_deg = 15.0 * (local_solar_time - 12.0);
        real hour_angle_H_rad = hour_angle_H_deg * GlobalConstants::PI / 180.0;
        real declination_rad = 23.44 * GlobalConstants::PI / 180.0 * sin(B_rad);
        real cos_theta_z = sin(phi_rad) * sin(declination_rad) + cos(phi_rad) * cos(declination_rad) * cos(hour_angle_H_rad);
        real theta_z_rad = acos(std::max(-1.0, std::min(1.0, cos_theta_z)));
        if (abs(sin(theta_z_rad)) < 1e-6) return {0, 0, cos_theta_z};
        real sin_azimuth = -sin(hour_angle_H_rad) * cos(declination_rad) / sin(theta_z_rad);
        real cos_azimuth = (sin(declination_rad) * cos(phi_rad) - cos(declination_rad) * sin(phi_rad) * cos(hour_angle_H_rad)) / sin(theta_z_rad);
        return {sin(theta_z_rad) * sin_azimuth, sin(theta_z_rad) * cos_azimuth, cos_theta_z };
    }

    // --- AtmosphereCore Implementation ---
    AtmosphereCore::State AtmosphereCore::query_state_at(const Vector3& position) const {
        State s;
        real effective_alt_m = std::min(position.z, 11000.0);
        s.T_air_C = sea_level_T_air_C + lapse_rate_C_per_m * effective_alt_m;
        s.p_air_kPa = sea_level_p_air_kPa * pow(1 + (lapse_rate_C_per_m / (sea_level_T_air_C + 273.15)) * effective_alt_m, -GlobalConstants::G / (lapse_rate_C_per_m * GlobalConstants::R_d));
        s.RH_percent = sea_level_RH_percent;
        real T_K = s.T_air_C + 273.15, p_Pa = s.p_air_kPa * 1000.0;
        real es_Pa = 0.61094 * exp((17.625 * s.T_air_C) / (s.T_air_C + 243.04)) * 1000;
        real e_Pa = (s.RH_percent / 100.0) * es_Pa;
        s.rho_air_kg_m3 = (p_Pa - e_Pa) / (GlobalConstants::R_d * T_K) + e_Pa / (GlobalConstants::R_v * T_K);
        return s;
    }

    // --- StochasticEventManager Implementation ---
    void StochasticEventManager::initialize(unsigned int seed) {
        rng.seed(seed);
        event_type_dist = std::discrete_distribution<>({0.5, 0.5}); // gust, rain_cell
    }
    void StochasticEventManager::add_event(Event new_event) {
        active_events.push_back(new_event);
    }
    void StochasticEventManager::update(real dt_s, real current_sim_time, const WindAndTurbulence& wind_module, const TrackGeometry& track_geo, const Vector3& grid_dims_m) {
        time_since_last_event += dt_s;
        if (uniform_dist(rng) < (1.0 - exp(-poisson_rate_events_per_second * time_since_last_event))) {
            time_since_last_event = 0.0;
            int event_type_idx = event_type_dist(rng);
            Event new_event;
            new_event.start_time = current_sim_time;
            
            // Spawn events near track sectors if available
            auto sectors = track_geo.get_sectors();
            if (!sectors.empty()) {
                int sector_idx = std::uniform_int_distribution<>(0, sectors.size() - 1)(rng);
                new_event.position = sectors[sector_idx].center;
                new_event.position.x += (uniform_dist(rng) - 0.5) * sectors[sector_idx].radius;
                new_event.position.y += (uniform_dist(rng) - 0.5) * sectors[sector_idx].radius;
            } else {
                 new_event.position = { (uniform_dist(rng) - 0.5) * grid_dims_m.x, (uniform_dist(rng) - 0.5) * grid_dims_m.y, 0.0 };
            }

            new_event.velocity = wind_module.get_mean_wind_at_ref_height();

            if (event_type_idx == 0) { // gust
                new_event.type = "gust";
                new_event.duration = 10.0 + uniform_dist(rng) * 50.0;
                new_event.params["amplitude_m_s"] = 5.0 + uniform_dist(rng) * 10.0;
                new_event.params["radius_m"] = 50.0 + uniform_dist(rng) * 100.0;
                active_events.push_back(new_event);
            } else if (event_type_idx == 1) { // rain_cell
                new_event.type = "rain_cell";
                new_event.duration = 300.0 + uniform_dist(rng) * 600.0;
                new_event.params["intensity_mm_hr"] = 5.0 + uniform_dist(rng) * 45.0;
                new_event.params["radius_m"] = 200.0 + uniform_dist(rng) * 800.0;
                active_events.push_back(new_event);
            }
        }
        for(auto& event : active_events) {
            event.position.x += event.velocity.x * dt_s; event.position.y += event.velocity.y * dt_s;
        }
        active_events.erase(std::remove_if(active_events.begin(), active_events.end(), 
            [current_sim_time](const Event& e) { return (current_sim_time > e.start_time + e.duration); }), active_events.end());
    }

    // --- CloudModel Implementation ---
    void CloudModel::initialize(unsigned int seed) { cloud_noise = std::make_unique<PerlinNoise>(seed); }
    void CloudModel::update(real dt_s, const Vector3& wind_advection) {
        advection_offset.x += wind_advection.x * dt_s;
        advection_offset.y += wind_advection.y * dt_s;
    }
    real CloudModel::get_cloud_cover_at(const Vector3& position, real sim_time) const {
        real scale = 0.0002;
        real noise_val = cloud_noise->noise((position.x + advection_offset.x) * scale, (position.y + advection_offset.y) * scale, sim_time * 0.01);
        return std::max(0.0, std::min(1.0, (noise_val + 0.2) / 1.2));
    }

    // --- SolarRadiation Implementation ---
    SolarRadiation::Irradiance SolarRadiation::query_irradiance(const Vector3& position, const Vector3& sun_vector, const AtmosphereCore::State& air_state, real local_cloud_cover, const TrackGeometry& track_geo) const {
        Irradiance result = {0.0, 0.0, 0.0, 90.0, false};
        real cos_zenith = sun_vector.z;

        if (cos_zenith <= 0) {
            result.zenith_angle_deg = 90.0;
            return result; 
        }

        for(const auto& caster : track_geo.get_shadow_casters()) {
             // Simplified AABB shadow check
            real shadow_length = 100; // Assume height of object casting shadow
            Vector3 shadow_end = {caster.min_corner.x - sun_vector.x * shadow_length, caster.min_corner.y - sun_vector.y * shadow_length, 0};
            if(position.x > caster.min_corner.x && position.x < caster.max_corner.x && position.y > caster.min_corner.y && position.y < caster.max_corner.y) {
                 result.is_in_shadow = true;
            }
        }
        
        real zenith_rad = acos(cos_zenith);
        result.zenith_angle_deg = zenith_rad * 180.0 / GlobalConstants::PI;
        real air_mass = 1.0 / (cos_zenith + 0.50572 * pow(96.07995 - result.zenith_angle_deg, -1.6364));
        real clear_sky_transmittance = pow(0.7, pow(air_mass, 0.678));
        real cloud_factor = 1.0 - 0.75 * pow(local_cloud_cover, 3.4);
        real transmittance = clear_sky_transmittance * cloud_factor;
        
        result.I_direct_W_m2 = result.is_in_shadow ? 0.0 : GlobalConstants::SOLAR_CONSTANT * transmittance * cos_zenith;
        result.I_diffuse_W_m2 = GlobalConstants::SOLAR_CONSTANT * (0.271 - 0.294 * clear_sky_transmittance) * cos_zenith * (1.0 + 0.5 * local_cloud_cover);
        
        real T_air_K = air_state.T_air_C + 273.15;
        real sky_emissivity = (0.72 + 0.005 * air_state.T_air_C) * (1.0 + 0.22 * pow(local_cloud_cover, 2));
        result.I_longwave_down_W_m2 = sky_emissivity * GlobalConstants::SIGMA_SB * pow(T_air_K, 4);
        return result;
    }

    // --- WindAndTurbulence Implementation ---
    void WindAndTurbulence::initialize(unsigned int seed) { turbulence_noise = std::make_unique<PerlinNoise>(seed); }
    WindAndTurbulence::WindState WindAndTurbulence::query_wind_at(const Vector3& position, real sim_time_s, const std::vector<StochasticEventManager::Event>& events) const {
        WindState state;
        real z = std::max(surface_roughness_z0, position.z);
        real speed_at_z = (friction_velocity_u_star / GlobalConstants::VON_KARMAN_K) * log(z / surface_roughness_z0);
        real ref_speed_mag = sqrt(ref_wind_speed.x*ref_wind_speed.x + ref_wind_speed.y*ref_wind_speed.y);
        state.wind_vector_m_s = ref_speed_mag > 0.01 ? Vector3{speed_at_z * (ref_wind_speed.x/ref_speed_mag), speed_at_z * (ref_wind_speed.y/ref_speed_mag), 0.0} : Vector3{0,0,0};
        
        state.turbulence_sigma_u = 2.5 * friction_velocity_u_star * (1.0 - 0.8 * (z / 500.0));
        real scale = 0.05;
        state.wind_vector_m_s.x += turbulence_noise->noise(position.x*scale, position.y*scale, sim_time_s*0.1) * state.turbulence_sigma_u;
        state.wind_vector_m_s.y += turbulence_noise->noise(position.y*scale+10.0, position.x*scale, sim_time_s*0.1) * state.turbulence_sigma_u;
        state.wind_vector_m_s.z += turbulence_noise->noise(position.z*scale, sim_time_s*0.1, position.x*scale+20.0) * state.turbulence_sigma_u * 0.5;

        for(const auto& event : events) {
            if(event.type == "gust") {
                real dist_sq = pow(position.x - event.position.x, 2) + pow(position.y - event.position.y, 2);
                real radius = event.params.at("radius_m");
                if (dist_sq < pow(radius, 2)) {
                    real dist_factor = 1.0 - sqrt(dist_sq) / radius;
                    real time_in_event = sim_time_s - event.start_time;
                    if (event.duration > 0 && time_in_event > 0) {
                        real time_factor = sin((time_in_event / event.duration) * GlobalConstants::PI);
                        real gust_magnitude = event.params.at("amplitude_m_s") * dist_factor * time_factor;
                        if (ref_speed_mag > 0.01) {
                            state.wind_vector_m_s.x += gust_magnitude * (ref_wind_speed.x / ref_speed_mag);
                            state.wind_vector_m_s.y += gust_magnitude * (ref_wind_speed.y / ref_speed_mag);
                        }
                    }
                }
            }
        }
        return state;
    }
    
    // --- PrecipitationAndHydrometeors Implementation ---
    PrecipitationAndHydrometeors::PrecipData PrecipitationAndHydrometeors::query_rain_at(const Vector3& position, const std::vector<StochasticEventManager::Event>& events) const {
        PrecipData data;
        for(const auto& event : events) {
            if(event.type == "rain_cell") {
                real dist_sq = pow(position.x - event.position.x, 2) + pow(position.y - event.position.y, 2);
                if (dist_sq < pow(event.params.at("radius_m"), 2)) {
                    data.rain_intensity_mm_hr += event.params.at("intensity_mm_hr");
                }
            }
        }
        if (data.rain_intensity_mm_hr > 0.1) {
            data.mean_drop_diameter_mm = 1.3 * pow(data.rain_intensity_mm_hr, 0.232);
            data.dsd_shape_param = 0.0;
        }
        return data;
    }

    // --- MoistureAndHydrology Implementation ---
    void MoistureAndHydrology::update(real dt_s, const PrecipitationAndHydrometeors::PrecipData& precip, const AtmosphereCore::State& air, 
                                      const SurfaceThermodynamics& surface, const WindAndTurbulence::WindState& wind) {
        real precip_m_s = (precip.rain_intensity_mm_hr / 1000.0) / 3600.0;
        surface_water_depth_mm += precip_m_s * dt_s * 1000.0;
        
        real T_surf_C = surface.get_surface_temp_C();
        real es_Pa = 0.61094 * exp((17.625 * T_surf_C) / (T_surf_C + 243.04)) * 1000;
        real ea_Pa = (air.RH_percent / 100.0) * 0.61094 * exp((17.625 * air.T_air_C) / (air.T_air_C + 243.04)) * 1000;
        real vpd_Pa = std::max(0.0, es_Pa - ea_Pa);
        real wind_speed = sqrt(pow(wind.wind_vector_m_s.x, 2) + pow(wind.wind_vector_m_s.y, 2));
        
        real evaporation_m_s = GlobalConstants::EVAPORATION_COEFF * (1.0 + 0.5 * wind_speed) * vpd_Pa;
        real actual_evap_mm = std::min(surface_water_depth_mm, evaporation_m_s * dt_s * 1000.0);
        surface_water_depth_mm -= actual_evap_mm;
        current_flux.evaporation_mm_s = actual_evap_mm > 0 ? (actual_evap_mm/1000.0) / dt_s : 0.0;
        
        real max_infil_rate_m_s = GlobalConstants::MAX_INFILTRATION_RATE_M_S;
        real infil_rate_m_s = max_infil_rate_m_s * (1.0 - soil_moisture_content);
        real actual_infil_mm = std::min(surface_water_depth_mm, infil_rate_m_s * dt_s * 1000.0);
        surface_water_depth_mm -= actual_infil_mm;
        soil_moisture_content = std::min(1.0, soil_moisture_content + actual_infil_mm * 0.005);
        
        current_flux.precipitation_mm_s = precip_m_s;
    }

    // --- SurfaceThermodynamics Implementation ---
    void SurfaceThermodynamics::update(real dt_s, const SolarRadiation::Irradiance& solar, const AtmosphereCore::State& air, const MoistureAndHydrology& hydro) {
        real T_surf_K = surface_temp_C + 273.15;
        real solar_absorbed = (solar.I_direct_W_m2 + solar.I_diffuse_W_m2) * (1.0 - albedo);
        real longwave_out = emissivity * GlobalConstants::SIGMA_SB * pow(T_surf_K, 4);
        real convective_out = GlobalConstants::CONVECTIVE_HEAT_COEFF * (surface_temp_C - air.T_air_C);
        real evaporative_cooling = hydro.get_fluxes().evaporation_mm_s * GlobalConstants::WATER_DENSITY * GlobalConstants::LATENT_HEAT_VAPORIZATION_LV;
        real net_flux_W_m2 = solar_absorbed + solar.I_longwave_down_W_m2 - longwave_out - convective_out - evaporative_cooling;
        surface_temp_C += (net_flux_W_m2 / heat_capacity_J_m2_K) * dt_s;
    }

    // --- AerosolsAndVisibility Implementation ---
    void AerosolsAndVisibility::update(real dt_s, const WindAndTurbulence::WindState& wind, const PrecipitationAndHydrometeors::PrecipData& precip) {
        real wind_speed_sq = pow(wind.wind_vector_m_s.x, 2) + pow(wind.wind_vector_m_s.y, 2);
        real resuspension_ug_s = 0.005 * wind_speed_sq;
        aerosol_conc_ug_m3 += resuspension_ug_s * dt_s;
        real washout_fraction_s = 0.0005 * precip.rain_intensity_mm_hr;
        aerosol_conc_ug_m3 *= (1.0 - washout_fraction_s * dt_s);
        aerosol_conc_ug_m3 *= (1.0 - 0.00002 * dt_s);
        aerosol_conc_ug_m3 = std::max(5.0, aerosol_conc_ug_m3);
    }
    real AerosolsAndVisibility::get_visibility_m(const Vector3& position, real rh_percent, const std::vector<StochasticEventManager::Event>& events) const {
        real b_ext_aerosol = 0.00003 * aerosol_conc_ug_m3;
        real rh_factor = 1.0;
        if (rh_percent > 90.0) {
            rh_factor = 1.0 + 4.0 * pow((rh_percent - 90.0) / 10.0, 2.0);
        }
        real visibility = 3.912 / std::max(1e-6, b_ext_aerosol * rh_factor);

        // Vehicle spray not implemented in this version
        return visibility;
    }

    // --- SurfaceGrid Implementation ---
    void SurfaceGrid::initialize(int x_dim, int y_dim, real cell_size_m, const Vector3& origin) {
        grid_x_dim = x_dim;
        grid_y_dim = y_dim;
        cell_size = cell_size_m;
        grid_origin = origin;

        grid.resize(y_dim, std::vector<Cell>(x_dim));
        for(int y = 0; y < y_dim; ++y) {
            for(int x = 0; x < x_dim; ++x) {
                grid[y][x].center_pos = {
                    grid_origin.x + (x + 0.5) * cell_size,
                    grid_origin.y + (y + 0.5) * cell_size,
                    0.0
                };
                // Create a slight depression in the middle for water to pool
                real dx = (x - x_dim/2.0);
                real dy = (y - y_dim/2.0);
                grid[y][x].altitude = 0.0001 * (dx*dx + dy*dy);
                grid[y][x].depression_depth_mm = 2.0;
            }
        }
    }
    int SurfaceGrid::get_x_idx(real world_x) const {
        int idx = static_cast<int>(floor((world_x - grid_origin.x) / cell_size));
        return std::max(0, std::min(grid_x_dim - 1, idx));
    }
    int SurfaceGrid::get_y_idx(real world_y) const {
        int idx = static_cast<int>(floor((world_y - grid_origin.y) / cell_size));
        return std::max(0, std::min(grid_y_dim - 1, idx));
    }

    SurfaceGrid::Cell& SurfaceGrid::get_cell_at(const Vector3& position) {
        return grid[get_y_idx(position.y)][get_x_idx(position.x)];
    }

    const SurfaceGrid::Cell& SurfaceGrid::get_cell_at(const Vector3& position) const {
         return grid[get_y_idx(position.y)][get_x_idx(position.x)];
    }


    void SurfaceGrid::update_all_physics(real dt_s, const Environment& env) {
        for(auto& row : grid) {
            for(auto& cell : row) {
                auto air_state = env.get_atmosphere_core().query_state_at(cell.center_pos);
                auto wind_state = env.get_wind_module().query_wind_at(cell.center_pos, env.get_time_manager().get_sim_time(), env.get_event_manager().get_active_events());
                auto precip_state = env.get_precipitation_module().query_rain_at(cell.center_pos, env.get_event_manager().get_active_events());
                auto solar_state = env.querySolar(cell.center_pos);
                
                cell.hydro_module.update(dt_s, precip_state, air_state, cell.thermo_module, wind_state);
                cell.thermo_module.update(dt_s, solar_state, air_state, cell.hydro_module);
            }
        }
    }

    void SurfaceGrid::update_water_flow(real dt_s) {
        if(grid.empty()) return;

        std::vector<std::vector<real>> net_flow_mm(grid_y_dim, std::vector<real>(grid_x_dim, 0.0));

        for(int y = 0; y < grid_y_dim; ++y) {
            for(int x = 0; x < grid_x_dim; ++x) {
                auto& current_cell = grid[y][x];
                real current_total_head = current_cell.altitude + current_cell.hydro_module.get_surface_water_depth_mm() / 1000.0;
                
                // Check neighbors
                int neighbors[4][2] = {{x, y - 1}, {x, y + 1}, {x - 1, y}, {x + 1, y}};
                for(auto& n : neighbors) {
                    int nx = n[0]; int ny = n[1];
                    if (nx >= 0 && nx < grid_x_dim && ny >= 0 && ny < grid_y_dim) {
                        auto& neighbor_cell = grid[ny][nx];
                        real neighbor_total_head = neighbor_cell.altitude + neighbor_cell.hydro_module.get_surface_water_depth_mm() / 1000.0;
                        real head_difference = current_total_head - neighbor_total_head;
                        
                        real flow_rate_m3_s = 0.1 * cell_size * head_difference; // Simplified flow
                        real flow_volume_m3 = flow_rate_m3_s * dt_s;
                        real flow_depth_change_m = flow_volume_m3 / (cell_size * cell_size);
                        
                        net_flow_mm[y][x] -= flow_depth_change_m * 1000.0;
                    }
                }
            }
        }

        for(int y = 0; y < grid_y_dim; ++y) {
            for(int x = 0; x < grid_x_dim; ++x) {
                grid[y][x].hydro_module.add_surface_water(net_flow_mm[y][x]);
            }
        }
    }

    // --- Environment (API)  ---
    void Environment::initialize(const std::map<std::string, std::any>& config, TrackGeometry track) {
        unsigned int seed = 1337;
        if(config.count("seed")) seed = std::any_cast<unsigned int>(config.at("seed"));
        
        if(config.count("latitude_deg")) latitude_deg = std::any_cast<real>(config.at("latitude_deg"));
        if(config.count("longitude_deg")) longitude_deg = std::any_cast<real>(config.at("longitude_deg"));
        if(config.count("start_time_iso")) time_manager.set_time(std::any_cast<std::string>(config.at("start_time_iso")));
        
        int grid_dim_x = 10, grid_dim_y = 10;
        real cell_size = 100.0;
        Vector3 grid_origin = {-500.0, -500.0, 0.0};
        if(config.count("grid_dim_x")) grid_dim_x = std::any_cast<int>(config.at("grid_dim_x"));
        if(config.count("grid_dim_y")) grid_dim_y = std::any_cast<int>(config.at("grid_dim_y"));
        if(config.count("cell_size_m")) cell_size = std::any_cast<real>(config.at("cell_size_m"));
        if(config.count("grid_origin_x")) grid_origin.x = std::any_cast<real>(config.at("grid_origin_x"));
        if(config.count("grid_origin_y")) grid_origin.y = std::any_cast<real>(config.at("grid_origin_y"));
        
        surface_grid.initialize(grid_dim_x, grid_dim_y, cell_size, grid_origin);
        stochastic_event_manager.initialize(seed);
        wind_and_turbulence.initialize(seed + 1);
        cloud_model.initialize(seed + 2);
        track_geometry = track;

        track_condition_rates.resize(grid_dim_y, std::vector<TrackConditionRate>(grid_dim_x));

        initialized = true;
        std::cout << "Environment initialized for Lat/Lon: " << latitude_deg << "/" << longitude_deg << " with a " << grid_dim_x << "x" << grid_dim_y << " grid." << std::endl;
    }
    void Environment::step(real dt_s, const std::vector<VehicleState>& vehicles) {
        if (!initialized) throw std::runtime_error("Environment not initialized.");
        if (dt_s <= 0) return;

        time_manager.advance(dt_s);
        
        auto mean_wind = wind_and_turbulence.get_mean_wind_at_ref_height();
        cloud_model.update(dt_s, mean_wind);
        stochastic_event_manager.update(dt_s, time_manager.get_sim_time(), wind_and_turbulence, track_geometry, 
            {surface_grid.get_x_dim() * surface_grid.get_cell_size(), surface_grid.get_y_dim() * surface_grid.get_cell_size(), 0});
        
        for(auto& row : track_condition_rates) {
            std::fill(row.begin(), row.end(), TrackConditionRate{});
        }
        
        surface_grid.update_all_physics(dt_s, *this);
        surface_grid.update_water_flow(dt_s);
        
        // Update track condition rates based on vehicles
        for(const auto& vehicle : vehicles) {
            auto& cell = surface_grid.get_cell_at(vehicle.position);
            auto& rates = track_condition_rates[surface_grid.get_y_idx(vehicle.position.y)][surface_grid.get_x_idx(vehicle.position.x)];
            
            real temp_factor = std::max(0.0, std::min(1.0, (cell.thermo_module.get_surface_temp_C() - 15.0) / 20.0));
            rates.rubber_add_rate += vehicle.slip_factor * temp_factor * 0.001 * dt_s;
            if (vehicle.slip_factor > 0.8) {
                rates.marbles_formation_rate += (vehicle.slip_factor - 0.8) * 0.0005 * dt_s;
            }
        }
    }
    Environment::AtmosphereQueryResult Environment::queryAtmosphere(const Vector3& position) const {
        auto core_state = atmosphere_core.query_state_at(position);
        auto wind_state = wind_and_turbulence.query_wind_at(position, time_manager.get_sim_time(), stochastic_event_manager.get_active_events());
        auto visibility = visibility_module.get_visibility_m(position, core_state.RH_percent, stochastic_event_manager.get_active_events());
        return {core_state.T_air_C, core_state.p_air_kPa, core_state.RH_percent, wind_state.wind_vector_m_s, visibility};
    }
    Environment::SolarQueryResult Environment::querySolar(const Vector3& position) const {
        Vector3 sun_vec = time_manager.compute_sun_vector(latitude_deg, longitude_deg);
        auto air_state = atmosphere_core.query_state_at(position);
        real cloud_cover = cloud_model.get_cloud_cover_at(position, time_manager.get_sim_time());
        return solar_radiation.query_irradiance(position, sun_vec, air_state, cloud_cover, track_geometry);
    }
    Environment::PrecipitationQueryResult Environment::queryPrecipitation(const Vector3& position) const {
        return precipitation_module.query_rain_at(position, stochastic_event_manager.get_active_events());
    }
    Environment::SurfaceStateResult Environment::querySurfaceState(const Vector3& position) const {
        const auto& cell = surface_grid.get_cell_at(position);
        return {cell.thermo_module.get_surface_temp_C(), cell.hydro_module.get_surface_water_depth_mm()};
    }
    Environment::SurfaceModifierResult Environment::querySurfaceModifiers(const Vector3& position, real vehicle_speed_ms) const {
        const auto& cell = surface_grid.get_cell_at(position);
        SurfaceModifierResult result;

        real temp_factor = std::max(0.0, std::min(1.0, (cell.thermo_module.get_surface_temp_C() - 15.0) / 20.0));
        result.base_grip_modifier = 0.9 + std::max(0.0, temp_factor * 0.1);
        
        real water_depth_m = cell.hydro_module.get_surface_water_depth_mm() / 1000.0;
        if (water_depth_m > 0.001) { // 1mm
            result.rolling_resistance_delta = std::min(100.0, 50.0 * water_depth_m);
            // Aquaplaning model: risk increases with speed^2 and water depth
            real aquaplaning_risk = (vehicle_speed_ms * vehicle_speed_ms) * (water_depth_m * 10.0);
            aquaplaning_risk = std::max(0.0, std::min(1.0, aquaplaning_risk));
            result.base_grip_modifier *= (1.0 - aquaplaning_risk);
        }
        
        result.base_grip_modifier = std::max(0.1, std::min(1.1, result.base_grip_modifier));
        return result;
    }

    const std::vector<std::vector<Environment::TrackConditionRate>>& Environment::get_track_condition_rates() const {
        return track_condition_rates;
    }
}


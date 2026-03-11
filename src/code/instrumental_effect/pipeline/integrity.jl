"""
    _integrity_push_check!(...)

Adds one integrity check record to the report lists.
"""
function _integrity_push_check!(checks::Vector{Dict{String,Any}}, critical_failures::Vector{String},
                                warnings::Vector{String}, name::String, passed::Bool;
                                critical::Bool=true, details::AbstractString="")
    level = critical ? "critical" : "warning"
    msg = passed ? "ok" : (isempty(details) ? "failed" : details)
    push!(checks, Dict{String,Any}(
        "name" => name,
        "level" => level,
        "passed" => passed,
        "message" => msg,
    ))
    if !passed
        if critical
            push!(critical_failures, name)
        else
            push!(warnings, name)
        end
    end
    return nothing
end

"""
    _write_integrity_report(...)

Writes a compact integrity report next to run outputs.
"""
function _write_integrity_report(path::AbstractString, report::Dict{String,Any})
    lines = String[
        "integrity_status=$(report["status"])",
        "critical_failures=$(report["critical_failures"])",
        "warnings=$(report["warnings"])",
        "checks_total=$(report["checks_total"])",
        "",
        "[checks]",
    ]

    for item in report["checks"]
        push!(lines, string(
            "- ", item["name"],
            " | level=", item["level"],
            " | passed=", item["passed"],
            " | message=", item["message"],
        ))
    end

    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end
    return path
end

"""
    _validate_filter_inputs(...)

Validates dimensions, NaN/Inf, units consistency and expected physical ranges.
Critical failures abort the run.
"""
function _validate_filter_inputs(cfg::InstrumentalConfig, Qdata, Udata, Pmax0, nuArray, PhiArray)
    checks = Dict{String,Any}[]
    critical_failures = String[]
    warnings = String[]

    try
        validate_instrumental_config!(cfg)
        _integrity_push_check!(checks, critical_failures, warnings, "config.basics", true; critical=true)
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "config.basics", false;
                               critical=true, details=sprint(showerror, err))
    end

    nν = length(nuArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.non_empty", nν > 0; critical=true)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.phi.non_empty", length(PhiArray) > 0; critical=true)

    nu_finite = !isempty(nuArray) && all(isfinite, nuArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.finite", nu_finite; critical=true)
    nu_positive = !isempty(nuArray) && all(>(0), nuArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.positive_hz", nu_positive; critical=true)
    nu_sorted = length(nuArray) <= 1 || all(diff(nuArray) .> 0)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.nu.strictly_increasing", nu_sorted; critical=true)

    phi_finite = !isempty(PhiArray) && all(isfinite, PhiArray)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.phi.finite", phi_finite; critical=true)
    phi_sorted = length(PhiArray) <= 1 || all(diff(PhiArray) .> 0)
    _integrity_push_check!(checks, critical_failures, warnings, "axes.phi.strictly_increasing", phi_sorted; critical=true)

    if !isempty(nuArray)
        expected_min_hz = cfg.νmin_MHz * 1e6
        expected_max_hz = cfg.νmax_MHz * 1e6
        expected_step_hz = cfg.Δν_MHz * 1e6

        first_ok = isapprox(first(nuArray), expected_min_hz; rtol=1e-10, atol=1e-6)
        step_ok = length(nuArray) <= 1 || all(isapprox.(diff(nuArray), expected_step_hz; rtol=1e-10, atol=1e-6))
        # Range syntax (νmin:Δν:νmax) may stop before νmax when Δν does not divide the span.
        last_not_above_max = last(nuArray) <= expected_max_hz + 1e-6

        nν = length(nuArray)
        expected_last_from_step = expected_min_hz + (nν - 1) * expected_step_hz
        last_consistent_with_step = isapprox(last(nuArray), expected_last_from_step; rtol=1e-10, atol=1e-6)

        units_ok = first_ok && step_ok && last_not_above_max && last_consistent_with_step
        units_details = units_ok ? "" : string(
            "nu-axis mismatch (first=", first(nuArray),
            ", last=", last(nuArray),
            ", step_expected=", expected_step_hz,
            ", max_cfg=", expected_max_hz, ")"
        )
        _integrity_push_check!(checks, critical_failures, warnings, "units.frequency.MHz_to_Hz", units_ok;
                               critical=true, details=units_details)
    else
        _integrity_push_check!(checks, critical_failures, warnings, "units.frequency.MHz_to_Hz", false;
                               critical=true, details="empty nuArray")
    end

    try
        q_layout = require_stokes_cube_layout(Qdata, "Qnu", cfg, nν)
        u_layout = require_stokes_cube_layout(Udata, "Unu", cfg, nν)
        layouts_match = q_layout == u_layout
        _integrity_push_check!(checks, critical_failures, warnings, "dims.stokes_layout", layouts_match;
                               critical=true, details=layouts_match ? "" : "Qnu/Unu layout mismatch")
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "dims.stokes_layout", false;
                               critical=true, details=sprint(showerror, err))
    end

    try
        require_channel_index(cfg.ichan, nν, "ichan")
        _integrity_push_check!(checks, critical_failures, warnings, "index.ichan.in_bounds", true; critical=true)
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "index.ichan.in_bounds", false;
                               critical=true, details=sprint(showerror, err))
    end

    try
        require_map_layout(Pmax0, "Pmax", cfg)
        _integrity_push_check!(checks, critical_failures, warnings, "dims.pmax_layout", true; critical=true)
    catch err
        _integrity_push_check!(checks, critical_failures, warnings, "dims.pmax_layout", false;
                               critical=true, details=sprint(showerror, err))
    end

    q_finite = all(isfinite, Qdata)
    u_finite = all(isfinite, Udata)
    pmax_finite = all(isfinite, Pmax0)
    _integrity_push_check!(checks, critical_failures, warnings, "nan_inf.Qnu", q_finite; critical=true)
    _integrity_push_check!(checks, critical_failures, warnings, "nan_inf.Unu", u_finite; critical=true)
    _integrity_push_check!(checks, critical_failures, warnings, "nan_inf.Pmax", pmax_finite; critical=true)

    pmax_non_negative = all(Pmax0 .>= 0)
    _integrity_push_check!(checks, critical_failures, warnings, "physical.Pmax.non_negative", pmax_non_negative;
                           critical=true, details="Pmax must be >= 0")

    freq_band_plausible = (1.0 <= cfg.νmin_MHz <= cfg.νmax_MHz <= 5.0e4)
    _integrity_push_check!(checks, critical_failures, warnings, "physical.frequency_band_plausible", freq_band_plausible;
                           critical=false, details="Expected ~[1, 50000] MHz")

    phi_span_plausible = !isempty(PhiArray) && (maximum(abs.(PhiArray)) <= 1.0e5)
    _integrity_push_check!(checks, critical_failures, warnings, "physical.phi_range_plausible", phi_span_plausible;
                           critical=false, details="Expected |phi| <= 1e5 rad m^-2")

    report = Dict{String,Any}(
        "status" => isempty(critical_failures) ? "pass" : "fail",
        "critical_failures" => length(critical_failures),
        "warnings" => length(warnings),
        "checks_total" => length(checks),
        "critical_failed_checks" => copy(critical_failures),
        "warning_checks" => copy(warnings),
        "checks" => checks,
    )

    report_path = _write_integrity_report(joinpath(cfg.base_out, "integrity_report.txt"), report)
    report["report_path"] = report_path

    if isempty(critical_failures)
        @info "Integrity checks passed" checks=length(checks) warnings=length(warnings) report=report_path
    else
        @error "Integrity critical checks failed; aborting run" failed=critical_failures report=report_path
        error("Critical integrity checks failed: $(join(critical_failures, ", "))")
    end

    return report
end

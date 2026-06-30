# Reporte EDA Falcon Challenge
## Archivos analizados
- `DataSetExport-Discharge.Best Available@08461300-Instantaneous-m^3 s-20260629185451.csv` → `Discharge.Best Available@08461300`, estación `08461300`, tipo `R_obs_release`, unidad `m^3/s`
- `DataSetExport-Lake Area.Best Available@08461200-Instantaneous-m^2-20260629185344.csv` → `Lake Area.Best Available@08461200`, estación `08461200`, tipo `lake_area`, unidad `m^2`
- `DataSetExport-Percentage.Conservation-Web-Telemetry@08461200-Instantaneous-%-20260630094142.csv` → `Percentage.Conservation-Web-Telemetry@08461200`, estación `08461200`, tipo `percentage_conservation`, unidad `%`
- `DataSetExport-Reservoir Elevation.Web-Daily-m@08461200-Instantaneous-m-20260629185508.csv` → `Reservoir Elevation.Web-Daily-m@08461200`, estación `08461200`, tipo `reservoir_elevation`, unidad `m`
- `DataSetExport-Total Storage.Web-Daily-tcm@08461200-Instantaneous-m^3-20260629185416.csv` → `Total Storage.Web-Daily-tcm@08461200`, estación `08461200`, tipo `S_obs_storage`, unidad `m^3`

## Contraste contra PDF
- **S_obs(t)** / `Total Storage.Web-Daily-tcm@08461200`: PRESENTE
- **DeltaS_obs(t)** / `Discharge.Total.Change-in-Storage@08461200`: FALTA OFICIAL; se puede DERIVAR desde Total Storage como aproximación
- **R_obs(t)** / `Discharge.Best Available@08461300`: PRESENTE
- **Elevation(t)** / `Reservoir Elevation.Web-Daily-m@08461200`: PRESENTE
- **Area(t)** / `Lake Area.Best Available@08461200`: PRESENTE
- **Evaporation(t)** / `Evaporation,accumltd.Daily Evaporation - mm@08461200`: FALTANTE
- **Conservation %(t)** / `Percentage.Conservation-Web-Telemetry@08461200`: PRESENTE
- **Smax** / `Falcon total conservation storage capacity`: FALTA OFICIAL; se puede usar max(S_obs) como proxy NO OFICIAL

## Diagnóstico
- Falta `Discharge.Total.Change-in-Storage@08461200`; `DeltaS_obs` fue derivado desde `Total Storage` si se construyó el dataset semanal.
- Falta `Smax` oficial del reservoir overview; si se calcula Smin con `max(S_obs)`, debe reportarse como proxy no oficial.
- Si la descarga está en m³/s, se convirtió a volumen semanal como `mean(flow)*604800`.

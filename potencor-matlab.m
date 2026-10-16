function potencor_main()
    % POTENCOR MATLAB Port - Function Based (R2018b Compatible)
    % Computes EEG power and correlations using built-in FFT and vectorized statistics.
    clc;
    disp('--- POTENCOR MATLAB ---');
    
    % 1. User Input 
    nomar1 = input('Nombre del primer archivo de nombres ? ', 's');
    nomar2 = input('Nombre del segundo archivo de nombres ? ', 's');
    mismoar = strcmp(nomar1, nomar2);
    
    % Read file lists (MATLAB 2018b compatible using textscan)
    fid1 = fopen(nomar1, 'r');
    if fid1 ~= -1
        tmp1 = textscan(fid1, '%s', 'Delimiter', '\n');
        archivos1 = tmp1{1};
        fclose(fid1);
        archivos1(cellfun('isempty', archivos1)) = []; % Remove empty lines
    else
        error('No se pudo abrir el archivo: %s', nomar1);
    end
    
    fid2 = fopen(nomar2, 'r');
    if fid2 ~= -1
        tmp2 = textscan(fid2, '%s', 'Delimiter', '\n');
        archivos2 = tmp2{1};
        fclose(fid2);
        archivos2(cellfun('isempty', archivos2)) = []; % Remove empty lines
    else
        error('No se pudo abrir el archivo: %s', nomar2);
    end
    
    numarch = min([length(archivos1), length(archivos2), 190]);
    
    tipo_in = input('Son archivos de [T]exto o de [E]nteros ? ', 's');
    formatoentero = startsWith(upper(tipo_in), 'E');
    
    simarc = input('Existe archivo de marcas (S/N) ? ', 's');
    marcas = startsWith(upper(simarc), 'S');
    
    numpts = input('Puntos por segmento de EEG [128,256,512 o 1024] ? ');
    puntos = numpts;
    if puntos == 128
        puntos = 256;
    end
    
    parapro = input('Segmentos para los promedios (0=todos) ? ');
    frecmu = input('Frecuencia de muestreo (en Hertz) ? ');
    if numpts == 128
        frecmu = frecmu * 2;
    end
    
    % 2. Band or Frequency setup
    tipo_analisis = input('Analisis por [F]recuencia o por [B]andas ? ', 's');
    siind = startsWith(upper(tipo_analisis), 'F');
    
    fac = (puntos * 128.0) / (1024.0 * frecmu);
    
    % Initialize variables
    frent = 0; Hz32 = 0; bands_matrix = [];
    
    if siind
        frent = 25;
        if puntos == 256 && frecmu == 512
            frent = 16;
        end
        Hz32 = fix(207 * fac);
    else
        % Use predefined standard EEG bands based on sampling factor
        dei = fix(12 * fac); des = fix(31 * fac);
        te1i = fix(32 * fac); te1s = fix(47 * fac);
        te2i = fix(48 * fac); te2s = fix(63 * fac);
        al1i = fix(64 * fac); al1s = fix(80 * fac);
        al2i = fix(81 * fac); al2s = fix(104 * fac);
        be1i = fix(105 * fac); be1s = fix(152 * fac);
        be2i = fix(153 * fac); be2s = fix(207 * fac);
        
        % Adjust limits to prevent overlap 
        if dei == 0, dei = 1; end
        if des == te1i, te1i = te1i + 1; end
        if te1s == te2i, te2i = te2i + 1; end
        if te2s == al1i, al1i = al1i + 1; end
        if al1s == al2i, al2i = al2i + 1; end
        if al2s == be1i, be1i = be1i + 1; end
        if be1s == be2i, be2i = be2i + 1; end
        Hz32 = be2s;
        
        % Package bands for passing to the export function
        bands_matrix = [
            dei, des;
            te1i, te1s;
            te2i, te2s;
            al1i, al1s;
            al2i, al2s;
            be1i, be1s;
            be2i, be2s
        ];
    end
    
    nomaux = input('Nombre comun de los archivos de salida (max. 6 car.) ? ', 's');
    nomcom = extractBefore([nomaux '.'], '.');
    
    % 3. Core Processing Loop
    for i = 1:numarch
        f1_name = archivos1{i};
        f2_name = archivos2{i};
        
        % Load signals using vectorized I/O
        [x_sig, y_sig] = read_data(f1_name, f2_name, formatoentero, numpts);
        
        % Truncate to match segments
        n_samples = min(length(x_sig), length(y_sig));
        n_segs = floor(n_samples / puntos);
        if parapro > 0
            n_segs = min(n_segs, parapro);
        end
        
        x_sig = x_sig(1:n_segs * puntos);
        y_sig = y_sig(1:n_segs * puntos);
        
        % Reshape into [puntos x n_segs] for vectorized FFT
        X_mat = reshape(x_sig, puntos, n_segs);
        Y_mat = reshape(y_sig, puntos, n_segs);
        
        % Built-in FFT 
        X_fft = fft(X_mat);
        Y_fft = fft(Y_mat);
        
        % Compute Power and Cross-Spectral Density
        facdiv = 4.0 / (puntos^2);
        Pxx = abs(X_fft(2:Hz32+1, :)).^2 * facdiv;
        Pyy = abs(Y_fft(2:Hz32+1, :)).^2 * facdiv;
        Pxy = real(X_fft(2:Hz32+1, :) .* conj(Y_fft(2:Hz32+1, :))) * facdiv;
        
        % Average across segments
        Sxx = mean(Pxx, 2);
        Syy = mean(Pyy, 2);
        Sxy = mean(Pxy, 2);
        
        % Compute correlations and extract band values
        if ~siind
            export_bands(Sxx, Syy, Sxy, f1_name, f2_name, n_segs, bands_matrix, nomcom);
        else
            export_freqs(Sxx, Syy, Sxy, f1_name, f2_name, n_segs, nomcom, frent);
        end
    end
    disp(' **** PROGRAMA TERMINADO ****');
end

function [x, y] = read_data(f1, f2, formatoentero, numpts)
    if formatoentero
        fid1 = fopen(f1, 'r', 'l');
        fid2 = fopen(f2, 'r', 'l');
        x = fread(fid1, inf, 'int16');
        y = fread(fid2, inf, 'int16');
        fclose(fid1); fclose(fid2);
    else
        x = load(f1);
        y = load(f2);
    end
    
    % Original interpolation step for 128 points
    if numpts == 128
        x = interp1(1:length(x), x, 1:0.5:length(x), 'linear')';
        y = interp1(1:length(y), y, 1:0.5:length(y), 'linear')';
    end
end

function export_bands(Sxx, Syy, Sxy, nom1, nom2, segprom, bands, nomcom)
    n_bands = size(bands, 1);
    pow_x = zeros(1, n_bands);
    pow_y = zeros(1, n_bands);
    cor_xy = zeros(1, n_bands);
    
    for b = 1:n_bands
        idx = bands(b, 1):bands(b, 2);
        pow_x(b) = sum(Sxx(idx));
        pow_y(b) = sum(Syy(idx));
        % Correlation logic
        sum_xy = sum(Sxy(idx));
        if pow_x(b) >= 1.0 && pow_y(b) >= 1.0
            cor_xy(b) = sum_xy / sqrt(pow_x(b) * pow_y(b));
        else
            cor_xy(b) = 0.0;
        end
    end
    
    % Write to .ABS, .REL, .TER, .TRA
    fid_abs = fopen([nomcom '.ABS'], 'a');
    fprintf(fid_abs, '%14.1f', pow_x);
    fprintf(fid_abs, '%14.1f  %s %d\n', sum(pow_x), nom1, segprom);
    fclose(fid_abs);
    
    % Fisher Z Transform for correlations 
    fisher_z = atanh(cor_xy);
    fisher_z(abs(cor_xy) >= 1.0) = 99.99; % Out of bounds handling
    
    fid_ter = fopen([nomcom '.TER'], 'a');
    fprintf(fid_ter, '%7.3f', cor_xy);
    fprintf(fid_ter, '%7.2f', fisher_z);
    fprintf(fid_ter, '  %s %s %d\n', nom1, nom2, segprom);
    fclose(fid_ter);
end

function export_freqs(Sxx, Syy, Sxy, nom1, nom2, segprom, nomcom, frent)
    fid_abf = fopen([nomcom '.ABF'], 'a');
    fprintf(fid_abf, '%14.1f', Sxx(1:frent));
    fprintf(fid_abf, '  %s %d\n', nom1, segprom);
    fclose(fid_abf);
end
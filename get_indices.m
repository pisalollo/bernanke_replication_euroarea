function [indices] = get_indices(full_list, target_tickers)
    % GET_INDICES Restituisce gli indici numerici dei ticker richiesti
    %
    % INPUTS:
    %   full_list:      Cell array o String array con TUTTI i nomi delle variabili del dataset
    %                   (es. dataset.Properties.VariableNames o raw_list)
    %   target_tickers: Cell array o String array con i ticker di cui vuoi l'indice
    %                   (es. ["GDP_EA", "HICPIN_EA"])
    %
    % OUTPUT:
    %   indices:        Vettore colonna con gli indici numerici corrispondenti.
    
    % 1. Standardizzazione in string array per robustezza
    full_list = string(full_list);
    target_tickers = string(target_tickers);
    
    % 2. Trova la corrispondenza
    % LIA è logico (Trovato/Non trovato), LOCB è l'indice in full_list
    [found_mask, loc] = ismember(target_tickers, full_list);
    
    % 3. Estrai solo gli indici validi (quelli trovati)
    % Nota: loc contiene 0 dove non ha trovato nulla, quindi filtriamo
    indices = loc(found_mask);
    
    % 4. Controllo errori (FONDAMENTALE in dataset grandi)
    if ~all(found_mask)
        missing_items = target_tickers(~found_mask);
        fprintf('\nATTENZIONE: I seguenti %d ticker non sono stati trovati:\n', length(missing_items));
        disp(missing_items);
        warning('Controlla se ci sono typos o suffissi diversi (es. _EA)');
    end
end
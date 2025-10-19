import React, { useState, useEffect } from 'react';
import { useLocation } from 'react-router-dom';
import MoleculeViewer from '../viewer/MoleculeViewer';
import './ViewerPage.css';

async function geminiFindPdbId(moleculeName) {
  const apiKey = process.env.REACT_APP_GEMINI_API_KEY;
  if (!apiKey) {
    console.warn("Gemini API key not found. Please set REACT_APP_GEMINI_API_KEY in .env file");
    return null;
  }

  try {
    const response = await fetch(`https://generativelanguage.googleapis.com/v1beta/models/gemini-pro:generateContent?key=${apiKey}`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        contents: [{
          parts: [{
            text: `Find the PDB ID for the protein/molecule: "${moleculeName}". Return only the PDB ID (4 characters) or "NOT_FOUND" if no match. Examples: insulin -> 1INS, hemoglobin -> 1HHO, myoglobin -> 1MBN.`
          }]
        }]
      })
    });

    if (!response.ok) {
      throw new Error(`API request failed: ${response.status}`);
    }

    const data = await response.json();
    const result = data.candidates?.[0]?.content?.parts?.[0]?.text?.trim();
    
    if (result && result !== "NOT_FOUND" && result.length === 4) {
      return result.toUpperCase();
    }
    
    return null;
  } catch (error) {
    console.error("Gemini API error:", error);
    return null;
  }
}

export default function ViewerPage() {
  const [pdbText, setPdbText] = useState('');
  const [fileName, setFileName] = useState('');
  const [searchName, setSearchName] = useState('');
  const [searchLoading, setSearchLoading] = useState(false);
  const [searchError, setSearchError] = useState('');
  const location = useLocation();

  useEffect(() => {
    if (location.state && location.state.cifText) {
      setPdbText(location.state.cifText);
      setFileName(location.state.fileName || 'example.cif');
      // Auto-load the structure
      setTimeout(() => {
        const content = location.state.cifText;
        const name = location.state.fileName || 'example.cif';
        function detectFormat(name, text) {
          if (name) {
            const ext = name.split('.').pop().toLowerCase();
            if (ext === 'cif' || ext === 'mmcif') return 'cif';
            if (ext === 'pdb' || ext === 'ent' || ext === 'txt') return 'pdb';
          }
          const t = (text || '').trim();
          if (!t) return 'pdb';
          if (t.startsWith('data_') || t.includes('_atom_site.')) return 'cif';
          if (t.split('\n')[0].startsWith('HEADER') || t.includes('ATOM') || t.includes('HETATM')) return 'pdb';
          return 'pdb';
        }
        const format = detectFormat(name, content);
        window.dispatchEvent(new CustomEvent('loadStructure', { detail: { content, name, format } }));
      }, 0);
    }
  }, [location.state]);

  function handleFile(e) {
    const f = e.target.files && e.target.files[0];
    if (!f) return;
    setFileName(f.name);
    const reader = new FileReader();
    reader.onload = (ev) => setPdbText(ev.target.result);
    reader.readAsText(f);
  }

  return (
    <div className="viewer-page">
      <aside className="controls">
        <h2>Load Structure</h2>
        {/* Search by molecule name using Gemini/RCSB */}
        <div className="search-section">
          <label htmlFor="mol-name-search">Find by Molecule Name (AI):</label>
          <input
            id="mol-name-search"
            type="text"
            value={searchName}
            onChange={e => setSearchName(e.target.value)}
            placeholder="e.g. hemoglobin, insulin"
            disabled={searchLoading}
          />
          <button
            className='search-btn'
            disabled={!searchName.trim() || searchLoading}
            onClick={async () => {
              setSearchLoading(true);
              setSearchError('');
              try {
                const pdbId = await geminiFindPdbId(searchName.trim());
                if (!pdbId) throw new Error('No matching structure found.');
                // Fetch mmCIF from RCSB
                const cifResp = await fetch(`https://files.rcsb.org/download/${pdbId}.cif`);
                if (!cifResp.ok) throw new Error('Failed to fetch mmCIF.');
                const cifText = await cifResp.text();
                setPdbText(cifText);
                setFileName(`${pdbId}.cif`);
                // Auto-load into viewer
                window.dispatchEvent(new CustomEvent('loadStructure', { detail: { content: cifText, name: `${pdbId}.cif`, format: 'cif' } }));
              } catch (err) {
                setSearchError(err.message || 'Search failed.');
              } finally {
                setSearchLoading(false);
              }
            }}
          >
            {searchLoading ? 'Searching...' : 'Find & Load'}
          </button>
          {searchError && <div className="search-error">{searchError}</div>}
        </div>
        <input type="file" accept=".pdb,.ent,.txt,.cif,.mmcif" onChange={handleFile} />
        <p>or paste structure file content below (PDB or CIF/mmCIF)</p>
        <textarea
          value={pdbText}
          onChange={(e) => setPdbText(e.target.value)}
          placeholder="Paste PDB or mmCIF content here"
        />
        <div className="buttons">
          <button onClick={() => {
            const content = pdbText;
            const name = fileName;
            // Get selected representation
            const repSelect = document.getElementById('representation');
            const representation = repSelect ? repSelect.value : 'cartoon';
            // basic format detection by filename or content
            function detectFormat(name, text) {
              if (name) {
                const ext = name.split('.').pop().toLowerCase();
                if (ext === 'cif' || ext === 'mmcif') return 'cif';
                if (ext === 'pdb' || ext === 'ent' || ext === 'txt') return 'pdb';
              }
              const t = (text || '').trim();
              if (!t) return 'pdb';
              if (t.startsWith('data_') || t.includes('_atom_site.')) return 'cif';
              if (t.split('\n')[0].startsWith('HEADER') || t.includes('ATOM') || t.includes('HETATM')) return 'pdb';
              return 'pdb';
            }

            const format = detectFormat(name, content);
            window.dispatchEvent(new CustomEvent('loadStructure', { detail: { content, name, format } }));
            window.dispatchEvent(new CustomEvent('changeRepresentation', { detail: representation }));
          }}>
            Load
          </button>
          <button onClick={() => {
            setPdbText('');
            setFileName('');
            window.dispatchEvent(new CustomEvent('clearPDB'));
          }}>
            Clear
          </button>
        </div>
        <div className="representation-dropdown">
          <label htmlFor="representation">Representation:</label>
          <select
            id="representation"
            onChange={(e) => {
              window.dispatchEvent(
                new CustomEvent('changeRepresentation', { detail: e.target.value })
              );
            }}
          >
            <option value="cartoon">Cartoon</option>
            <option value="stick">Stick</option>
            <option value="ballstick">Ball and Stick</option>
            <option value="surface">Surface</option>
            <option value="sphere">Sphere (CPK)</option>
            <option value="line">Line</option>
            <option value="cross">Cross</option>
          </select>
        </div>
      </aside>

      <section className="viewer-area">
        <MoleculeViewer />
      </section>
    </div>
  );
}

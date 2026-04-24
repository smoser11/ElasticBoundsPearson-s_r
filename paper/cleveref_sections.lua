-- cleveref_sections.lua
-- For LaTeX output: replace Quarto's "Section~\ref{sec-...}" with "\Cref{sec-...}"
-- so that cleveref's nameinlink option produces "Section X.X" as a single hyperlink.

function RawInline(el)
  if el.format ~= 'latex' then return el end
  local new_text = el.text:gsub('Section~\\ref%{(sec%-[^}]+)%}', '\\Cref{%1}')
  if new_text ~= el.text then
    return pandoc.RawInline('latex', new_text)
  end
  return el
end

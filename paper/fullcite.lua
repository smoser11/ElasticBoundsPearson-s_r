local refs = {}

function Div(div)
  -- Store references from bibliography
  if div.identifier then
    local ref_id = div.identifier:match('^ref%-(.*)$')
    if ref_id then
      refs[ref_id] = div.content
    end
  end
end

function Para(para)
  -- Look for "!fullcite bibkey"
  if #para.content > 0 and para.content[1].t == "Str" then
    local first_word = para.content[1].text
    if first_word:sub(1,9) == "!fullcite" then
      -- Get bibkey (assumes "!fullcite bibkey" in paragraph)
      local bibkey = nil
      -- Handle space or no space after "!fullcite"
      if #para.content > 1 and para.content[2].t == "Space" and #para.content > 2 and para.content[3].t == "Str" then
        bibkey = para.content[3].text
      elseif #para.content > 1 and para.content[2].t == "Str" then
        bibkey = para.content[2].text
      else
        bibkey = first_word:match("!fullcite%s+(.+)")
      end
      if bibkey and refs[bibkey] then
        return pandoc.Para(refs[bibkey])
      end
    end
  end
end

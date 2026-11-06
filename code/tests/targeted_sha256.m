function digest=targeted_sha256(path)
%TARGETED_SHA256 对文件原始字节计算SHA-256，用于冻结协议与断点续跑校验。
fid=fopen(path,'rb'); assert(fid>=0,'无法读取文件：%s',path);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8'); md.update(bytes);
end
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
